import SedmDb
import numpy as np
from astropy.io import fits

db = SedmDb.SedmDB()


def create_atomic_requests():
    """
    Finds requests without atomicrequets and creates their atomicrequests
    Returns:
        (0, "...") containing a list of the requests and which, if any, failed
    """
    # get the ids of the requests with no atomicrequest
    request_ids = db.execute_sql('SELECT id FROM request WHERE NOT EXISTS '
                                 '(SELECT id FROM atomicrequest WHERE request.id = atomicrequest.request_id')
    # TODO: make this only create for requests with status 'PENDING'/'ACTIVE'?
    failed = []
    for request in request_ids:
        success = create_request_atomic_requests(request)
        if success[0] == -1:
            failed.append(request)
    for request in failed:
        request_ids.remove(request)
    if failed:
        return (0, "Added atomicrequests for requests %s, failed to add for request  %s" % (request_ids, failed))
    else:
        return (0, "Added atomicrequests for requests %s" % (request_ids,))


def create_request_atomic_requests(request_id):
    """
    create atomicrequest entries for a given request
    Args:
        request_id: int
            id of the request that needs atomicrequests
    Returns:
        (-1, "ERROR: ...") if there was an issue
        (0, "Added (#) atomic requests for request (#)") if it was successful
    """
    status = db.execute_sql("SELECT status FROM request WHERE id=%s" % (request_id,))
    if not status:
        return (-1, "ERROR: request does not exist!")
    elif status == 'CANCELED':
        return (-1, "ERROR: request has been canceled!")
    elif status == 'EXPIRED':
        return (-1, "ERROR: request has expired!")

    if db.execute_sql("SELECT id FROM atomicrequest WHERE request_id='%s';" % (request_id,)):
        return (-1, "ERROR: atomicrequests already exist for that request!")
    request = db.execute_sql("SELECT object_id, exptime, priority, inidate, enddate,"
                             "cadence, phasesamples, sampletolerance, filters, nexposures, ordering "
                             "FROM request WHERE id='%s';" % (request_id,))[0]
    # TODO: implement cadence/phasesamples/sampletolerance (I have no idea how they interact with nexposures)
    pardic = {'object_id': request[0], 'priority': request[2], 'inidate': request[3],
              'enddate': request[4], 'request_id': request_id}
    obs_order = []
    if request[10]:
        for num_fil in request[10]:
            for n in range(int(num_fil[0])):  # the number should be single digit
                obs_order.append(num_fil[1:])
    elif request[9]:
        for filter_idx in range(len(request[8])):
            for n in range(request[9][filter_idx]):
                obs_order.append(request[8][filter_idx])
    else:
        return (-1, "ERROR: request contains neither nexposures nor ordering!")  # add_request should prevevnt this
    ifus = np.where(np.array(obs_order) == 'ifu')[0]
    if len(ifus) == 2:  # TODO: rewrite, have another way of indicating a/b
        obs_order[ifus[0]] = 'ifu_a'
        obs_order[ifus[1]] = 'ifu_b'
    if any([(filt not in ['u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b']) for filt in obs_order]):
        return (-1, "ERROR: either filters or ordering has an invalid entry!")
    for n, filter_des in enumerate(obs_order):
        pardic['filter'] = filter_des
        pardic['order_id'] = n
        # TODO: do exptime modifications here per filter and ab vs single exposure
        if 'ifu' in filter_des:
            pardic['exptime'] = request[1][0]  # TODO: make sure the sql returns the proper format
        else:
            pardic['exptime'] = request[1][1]
        add_return = db.add_atomic_request(pardic)
        if add_return[0] == -1:
            return (-1, "ERROR: adding atomicrequest (# %s, filter %s) failed." % (n + 1, filter_des))
    return (0, "Added %s atomic requests for request %s" % (len(obs_order), request_id))


def add_observation_fitsfile(fitsfile):
    """
    Adds an observation from a fitsfile, sets the atomicobservation to 'OBSERVED'. Checks if all the atomicobservations
    associated with its request are observed, if so it sets the request's status to 'COMPLETED'.
    Args:
        fitsfile: str (fits file path)
    """
    hdulist = fits.open(fitsfile)
    header = hdulist[0].header
    header_dict = {'mjd': header['JD']-2400000.5, 'airmass': header['AIRMASS'], 'exptime': header['EXPTIME'],
                   'fitsfile': fitsfile, 'lst': header['LST'], 'ra': ra_to_decimal(header['RA']),
                   'dec': dec_to_decimal(header['DEC']), 'tel_ra': header['TEL_RA'], 'tel_dec': header['TEL_DEC'],
                   'tel_az': header['TEL_AZ'], 'tel_el': header['TEL_EL'], 'tel_pa': header['TEL_PA'],
                   'ra_off': header['RA_OFF'], 'dec_off': header['DEC_OFF'], 'camera': header['CAM_NAME']}
    # header_dict['atomicrequest_id'] = header['ATOM_ID']
    # TODO: add ATOM_ID to the header (above)
    # TODO: add 'imtype', generate from fitsfile name?

    tel_stats = {'date': header['OBSDATE'], 'dome_status': header['DOMEST'], 'in_temp': header['IN_AIR'],
                 'in_humidity': header['IN_HUM'], 'in_dew': header['IN_DEW'], 'out_temp': header['OUT_AIR'],
                 'out_humidity': header['OUT_HUM'], 'out_dew': header['OUT_HUM'], 'wind_dir': header['WIND_DIR'],
                 'wsp_cur': header['WSP_CUR'], 'wsp_avg': header['WSP_AVG'], 'mir_temp': header['MIR_TEMP'],
                 'top_air': header['TOP_AIR'], 'pri_temp': header['PRI_TEMP'], 'sec_temp': header['SEC_TEMP'],
                 'flo_temp': header['FLO_TEMP'], 'bot_temp': header['BOT_TEMP'], 'mid_temp': header['MID_TEMP'],
                 'top_temp': header['TOP_TEMP']}
    obs_added = db.add_observation_fits(header_dict, tel_stats)
    # check if all of the request's atomicrequests are 'OBSERVED', if so set the request's status to 'COMPLETED'
    # TODO: must test with connected request/atomicrequest/fitsfile
    if obs_added[0] == 0:
        print db.update_atomic_request({'id': header_dict['atomicrequest_id'], 'status': 'OBSERVED'})  # better way than print?
        request_id = db.execute_sql("SELECT request_id FROM atomicrequest WHERE id='%s'"
                                    % (header_dict['atomicrequest_id'],))
        atomic_requests = db.get_request_atomicrequests(request_id)
        if all(request[5] == 'OBSERVED' for request in atomic_requests):
            db.update_request({'id': request_id, 'status': 'COMPLETED'})
    return obs_added


def ra_to_decimal(ra):
    hms = ra.split(':')
    return float(hms[0])*15+float(hms[1])/4+float(hms[2])/240


def dec_to_decimal(dec):
    dms = dec.split(':')
    return float(dms[0])+float(dms[1])/60+float(dms[2])/3600
