import os
import json
import glob
import requests
import commands
import datetime


# Path constants
add_target_url = 'http://nera.palomar.caltech.edu/cgi-bin/' \
                 'telescopes/p60/sedm/exFollowup.cgi'

pharos_spec_dir = '/scr2/sedmdrp/redux/'
pharos_phot_dir = '/scr2/sedm/phot/'

growth_base_url = 'http://skipper.caltech.edu:8080/cgi-bin/growth/'
growth_inst_url = growth_base_url + 'update_followup_config.cgi'
growth_stat_url = growth_base_url + 'update_followup_status.cgi'
growth_spec_url = growth_base_url + 'add_spec_auto.cgi'
growth_phot_url = growth_base_url + 'edit_phot_auto.cgi'
growth_view_source_url = growth_base_url + 'view_source.cgi?'

default_id = 65
user, pwd = open('/home/sedm/.growth_creds.txt', 'r').readlines()[0].split()


def write_json_file(pydict, output_file):
    """
    Write the python dictionary to a json file
    :param pydict: 
    :param output_file: 
    :return: json file path
    """

    jsonFile = open(output_file, 'w')
    jsonFile.write(json.dumps(pydict))
    jsonFile.close()

    return output_file


def timestamp():
    """
    UTC timestamp.  Use this when saving files
    :return: 
    """
    return datetime.datetime.utcnow().strftime("%Y%m%d_%H_%M_%S")


def send_instrument_configuration(instrument_id="",
                                  options_programs=None,
                                  post_url="", save=False):
    """
    Update the instrument configuration.  Updating this file will change the 
    available options shown on the growth marshal page
    
    :param instrument_id: 
    :param options_programs: 
    :param post_url: 
    :return: 
    """

    if not instrument_id:
        instrument_id = default_id

    if not post_url:
        post_url = add_target_url

    # 1. Create the main json dictionary file
    instConfig = {
        'instrument_id': instrument_id,
        'post_url': post_url,
        'options': []
    }

    # 2. Check the options dictionary for values if none use default for SEDm
    if options_programs:
        instConfig['options'] = options_programs
    else:
        fourshot = {'name': 'Followup',
                    'type': 'select',
                    'value': 'Four Shot (r,g,i,u)'}

        three_shot = {'name': 'Followup',
                      'type': 'select',
                      'value': 'Three Shot (r,g,i)'}

        ifu = {'name': 'Followup',
               'type': 'select',
               'value': 'IFU'}

        ifu_fourshot = {'name': 'Followup',
                        'type': 'select',
                        'value': 'Fourshot + IFU'}

        check1 = {"name": "Filters", "type": "check", "value": "r"}

        check2 = {"name": "Filters", "type": "check", "value": "g"}

        check3 = {"name": "Filters", "type": "check", "value": "i"}

        check4 = {"name": "Filters", "type": "check", "value": "ifu"}

        instConfig['options'] = [fourshot, three_shot, ifu, ifu_fourshot,
                                 check1, check2, check3, check4]

    # 3. Create a json file with the request and read it in to memory
    jsonFile = open(write_json_file(instConfig, 'config.txt'), 'r')

    # 4. Send json file request to growth marshal
    ret = requests.post(growth_inst_url,
                        auth=(user, pwd),
                        files={'jsonfile': jsonFile})

    # 5. Close the file and save it if needed
    jsonFile.close()

    if not save:
        os.remove('config.txt')

    # 6. Return the request response for user to determine if it was a success
    return ret


def update_request(status, request_id, instrument_id='',
                   output_dir='targets', filename='', save=True):
    """
    Function to update the status of any request as long as it has
    not been deleted. The new status will show up on the status section
    of the request on the growth marshal.
    
    :param filename: 
    :param status: 
    :param instrument_id: 
    :param request_id: 
    :param output_dir: 
    :param save: 
    :return: 
    """

    # 1. Make sure that we have the required fields
    if not request_id:
        print "Can't update a request without the request id"
        return False

    if not instrument_id:
        instrument_id = default_id

    if not filename:
        filename = "%s_%s.txt" % (str(request_id), timestamp())

    output_file = os.path.join(output_dir, filename)

    # 2. Create the new status dictionary
    statusConfig = {'instrument_id': instrument_id,
                    'request_id': request_id,
                    'new_status': status}

    # 3. Write and read in the json file to memory
    jsonFile = open(write_json_file(statusConfig, output_file), 'r')

    # 4. Send the request, close the file, and save if needed
    ret = requests.post(growth_stat_url, auth=(user, pwd),
                        files={'jsonfile': jsonFile})
    jsonFile.close()

    if not save:
        os.remove(output_file)

    # 5. Print the request response for the user
    print ret

    return True


def get_keywords_from_file(inputfile, keywords, sep=':'):
    """
    Get keywords from file.  It is dependent on files having a specific format
    where the keyword is on the left and the value on the right by some common
    seperator
    
    :param keywords: 
    :return: 
    """
    return_dict = {}
    for k, v in keywords.iteritems():
        out = commands.getoutput('grep %s %s' % (v, inputfile))
        if k.upper() == 'EXPTIME':
            outstr = out.split(sep,1)[-1]
            print outstr
            return_dict[k] = int(outstr)
        else:
            return_dict[k] = out.split(sep,1)[-1]
    return return_dict


def upload_phot(phot_file, instrument_id=65, request_id=''):
    """
    
    :param phot_file: 
    :param instrument_id: 
    :param request_id: 
    :return: 
    """

    with open(phot_file, 'r') as photometryFile:
        photometry = photometryFile.read()

    photometry = photometry.split('\n')
    columnNames = photometry[0].split(',')
    photometry = photometry[1:-1]

    photometryList = []
    for entry in photometry:
        newDict = {}
        photometryPoint = entry.split(',')
        for index, column in enumerate(columnNames):
            data = photometryPoint[index]
            if '"' in data:
                data = data.replace('"', '')
            elif data == 't':
                data = True
            elif data == 'f':
                data = False
            else:
                data = float(data)
            if data == 'None':
                data = None
            newDict[column] = data
        photometryList.append(newDict)

        submissionDict = {}
        submissionDict['photometry_list'] = photometryList
        submissionDict['instrument_id'] = instrument_id
        submissionDict['request_id'] = request_id

        jsonFile = open('photometryExample.txt', 'w')
        jsonFile.write(json.dumps(submissionDict))
        jsonFile.close()

        jsonFile = open('photometryExample.txt', 'r')
        ret = requests.post(growth_phot_url, auth=(user, pwd),
                            files={'jsonfile': jsonFile})
        jsonFile.close()
        return ret


def upload_spectra(spec_file, fill_by_file=False, instrument_id=65,
                   request_id='', exptime=3600, get_id_from_db=True,
                   observer='SEDmRobot', reducedby="Neill", obsdate="",
                   output_dir='targets/', format_type='ascii', save=True,
                   quality=1, check_quality=True, min_quality=2):
    """
    Add spectra to the growth marshal.  If the fill_by_file is selected then
    most of the keywords will be filled from the spectra file itself.  If
    the request has been canceled for some reason then it will not be possible
    to update the request
    
    :param fill_by_file: 
    :param request_id: 
    :param exptime: 
    :param observer: 
    :param reducedby: 
    :param obsdate: 
    :param output_dir: 
    :param format_type: 
    :param save: 
    :param instrument_id: 
    :param spec_file: 
    :return: 
    """
    if not request_id:
        print "Can't update without a request id"
        return False

    basefile = os.path.basename(spec_file)
    output = os.path.join(output_dir, basefile)
    output_file = output.replace('.txt', '.json')

    # 1. Create the mandatory keyword dictionary
    if not fill_by_file:
        submission_dict = {'exptime': exptime,
                           'obsdate': obsdate,
                           'reducedby': reducedby}
    else:
        keywords_dict = {'reducedby': 'REDUCER',
                         'obsdate': 'OBSUTC',
                         'exptime': 'EXPTIME',
                         'quality': 'QUALITY'}
        
        submission_dict = get_keywords_from_file(spec_file, keywords_dict)
        quality = int(submission_dict['quality'])
        del submission_dict['quality']

    print type(format_type), type(request_id), type(instrument_id)
    submission_dict.update({'format': format_type.rstrip().lstrip(),
                           'instrument_id': instrument_id,
                            'request_id': request_id,
                            'observer': observer.rstrip().lstrip()})
    
    print type(submission_dict['instrument_id'])

    # 1a. [Optional] Check the quality if it is smaller or equal to min quality
    # then upload spectra
    
    if check_quality and quality > min_quality:
        print "Spectra quality does not pass"
        return False
    
    # 2. Open the configuration and spec file for transmission
    jsonFile = open(write_json_file(submission_dict, output_file), 'r')
    upfile = open(spec_file, 'r')

    # 3. Send the request
    ret = requests.post(growth_spec_url, auth=(user, pwd),
                        files={'jsonfile': jsonFile, 'upfile': upfile})
    print ret

    # 4. Close files and send request response
    upfile.close()
    jsonFile.close()
    #if not save:
    #    os.remove(output_file)

    return True


def read_request(request_file):
    """
    Read in request file to python dictionary
    :param request_file: 
    :return: 
    """
    return json.load(open(request_file, 'r'))


def update_target_by_object(objname, add_spectra=False, spectra_file='',
                            add_status=False, status='Completed',
                            pull_requests=False,
                            add_phot=False, phot_file='', search_db=False,
                            target_dir='requests/', target_base_name='request'):
    """
    Go through the request and find the one that matches the objname
    :param objname: 
    :return: 
    """
    status_ret = False
    spec_ret = False
    phot_ret = False
    # 1. Start by looking at all files in the target directory
    # this is currently the directory that holds all incoming request
    # dictionary from growth.
    if pull_requests:
        print "Gathering all request files"
        import growth_watcher
        growth_watcher.pull_request_from_remote()
    match_list = []
    files = glob.glob('%s%s*' % (target_dir, target_base_name))

    # TODO: Change this to either look at the directory or to the SEDm database
    #if search_db:
    #    import SedmDb
    #    sedmdb = SedmDb.SedmDB("sedmdb", "pharos.caltech.edu")

    # 1a. Search for target in the database

    for i in files:
        targ = read_request(i)
        # Add any request file that contains the objname.  There may be more
        # than one if an update to the original request has been sent. Ignore
        # any request with the status delete as we can not update those
        if (targ['sourcename'].lower() == objname.lower() and
                    targ['status'] != 'delete'):

            match_list.append(targ)
         
    if len(match_list) == 1:
        print match_list
        target = match_list[0]
        print "Uploading target %s files" % objname
        if add_spectra:
            print target['requestid']
            spec_ret = upload_spectra(spectra_file, fill_by_file=True,
                                      request_id=target['requestid'])
        if add_phot:
            phot_ret = upload_phot(phot_file, request_id=target['requestid'])

        if add_status:
            status_ret = update_request(status, request_id=target['requestid'])

    elif len(match_list) == 0:
        print "Could not match name with any request file"

    else:
        print 'Multiple matches have been made for target: %s' % objname
        request_id_list = []
        for j in match_list:
            request_id_list.append(j['requestid'])
        
        # If all the request id matches then there is no problem and we can 
        # send data
        if all(x == request_id_list[0] for x in request_id_list):
            
            target = match_list[0]

        else:
            # TODO: Handle case in which an update has been sent
            # For now we just use highest value
            request_id = sorted(request_id_list)[-1]
            print request_id_list
            for j in match_list:
                if j['requestid'] == request_id:
                    target = j
                    break

        print "Updating target %s" % objname
        if add_spectra:
            print target['requestid']
            spec_ret = upload_spectra(spectra_file, fill_by_file=True,
                                      request_id=target['requestid'])
        if add_phot:
            phot_ret = upload_phot(phot_file, request_id=target['requestid'])

        if add_status:
            status_ret = update_request(status, request_id=target['requestid'])

    return_link = growth_view_source_url + "name=%s" % target['sourcename']

    return return_link, spec_ret, status_ret, phot_ret

          
def parse_ztf_by_dir(target_dir):
    """Given a target directory get all files that have ztf or ZTF as base 
       name"""

    if target_dir[-1] != '/':
        target_dir += '/'

    files = glob.glob('%sZTF*.txt' % target_dir)
    files += glob.glob('%sztf*.txt' % target_dir)

    out = open(target_dir + "report_ztf.txt", "w")
    out.write("ZTF Upload report for %s generated on %s\n\n" %
              (target_dir.split('/')[-2],
               datetime.datetime.now().strftime("%c")))
    pr = True
    for fi in files:
        objname = os.path.basename(fi).split('_')[0]
        r, spec, stat, phot = update_target_by_object(objname,
                                                      add_status=True,
                                                      status='Completed',
                                                      add_spectra=True,
                                                      spectra_file=fi,
                                                      pull_requests=pr)
        # Only need to pull requests the first time
        pr = False
        # print Object name
        out.write("%s: " % objname)
        # Was a spectrum uploaded?
        if spec:
            out.write("OK ")
        else:
            out.write("NO ")
        # Was status updated?
        if stat:
            out.write("OK ")
        else:
            out.write("NO ")
        print r
        out.write("%s\n" % r)

    # Close log file
    out.close()


if __name__ == '__main__':
    import os
    import datetime

    reddir = '/scr2/sedmdrp/redux/'
    utc = datetime.datetime.utcnow().strftime("%Y%m%d") 
    srcdir = reddir + utc + '/'
    if not os.path.exists(srcdir):
        print("Dir not found: %s" % srcdir)
    else:
        print("Uploading from %s" % srcdir)

        # Check requests, targets dirs
        reqdir = srcdir + 'requests'
        trgdir = srcdir + 'targets'
        if not os.path.exists(reqdir):
            os.mkdir(reqdir)
        if not os.path.exists(trgdir):
            os.mkdir(trgdir)
      
        parse_ztf_by_dir(srcdir)
        # Test to see if I can update from a given directory
