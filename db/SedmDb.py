import sqlalchemy.pool as pool
from sqlalchemy import exc, create_engine
import psycopg2
import numpy as np
import subprocess
import warnings
from astropy.io import fits


# Singleton/SingletonPattern.py

# noinspection PyArgumentList,SqlResolve,SqlNoDataSourceInspection,PyRedundantParentheses
class SedmDB:
    class __SedmDB:
        def __init__(self):
            """
            Creates the instance of db connections.
            Needs the username as a parameter.
            The password for SedmDB must be stored in ~/.pgpass
            """
            cmd = "cat ~/.pgpass | grep sedmdbtest | grep -v '#' | awk -F ':' '{print $4}' | head -n 1"

            p = subprocess.Popen(cmd, stdout=subprocess.PIPE, shell=True)
            self.user_sedmdb = p.stdout.read().replace('\n', '')

            self.pool_sedmdb = pool.QueuePool(self.__getSedmDBConn__, max_overflow=10, pool_size=2, recycle=True)

        def __str__(self):
            return repr(self)

        def __getSedmDBConn__(self):
            """
            Creates the connection to SedmDB.
            """
            sedmdbcon = psycopg2.connect(host="localhost", port="5432", dbname="sedmdbtest",
                                         user=self.user_sedmdb)
            return sedmdbcon

    instance = None

    def __init__(self):
        """
        Makes sure only one instance is created.
        """
        if not SedmDB.instance:
            SedmDB.instance = SedmDB.__SedmDB()

    def __getattr__(self, name):
        return getattr(self.instance, name)

    def execute_sql(self, sql):
        """
        Runs the SedmDB sql query in a safe way through the DBManager.

        Returns the object with the results.
        """
        conn = self.pool_sedmdb.connect()
        cursor = conn.cursor()
        try:
            cursor.execute(sql)
        except exc.DBAPIError, e:
            # an exception is raised, Connection is invalidated.
            if e.connection_invalidated:
                print "Connection was invalidated!"
        if 'SELECT' in sql:
            print cursor
            obj = cursor.fetchall()
            return obj
        else:
            cursor.execute("commit;")
            return []

    def get_conn_sedmDB(self):
        """
        Runs the WSDB sql query in a safe way through the DBManager.

        Returns the object with the results.
        """
        conn = self.pool_sedmdb.connect()

        return conn

    def add_user(self, pardic):
        """
        Adds a new user. Checks for duplicates in name. If user exists:
          (-1, "ERROR: User exists!")

        """
        # TODO: is name or e-mail the best designator to test for?
        # check if already contained
        users = self.execute_sql('SELECT * FROM users')
        usernames = [user[1] for user in users]
        # TODO: switch from group_id to designator?
        if not pardic['username'] in usernames:
            self.execute_sql("INSERT INTO users (username, name, email) Values ('%s', '%s', '%s');"
                             % (pardic['username'], pardic['name'], pardic['email']))
            user_id = self.execute_sql("SELECT id FROM users WHERE username='%s'" % (pardic['username'],))[0][0]
            print user_id, self.add_to_group(user_id, pardic['group_id'])
            return (0, "User added")
        else:
            return (-1, "ERROR: User with that username exists!")
        # TODO: test this

    def remove_user(self, pardic):
        """
        Removes the user. If user does not exist:
          (-1, "ERROR: User does not exist!")

        """
        # TODO: remove id option?
        if 'id' in pardic.keys():
            if pardic['id'] in [x[0] for x in self.execute_sql('SELECT id FROM users;')]:
                self.execute_sql("DELETE FROM usergroups WHERE id='%s'" % (user_id,))
                self.execute_sql("DELETE FROM users WHERE id='%s';" % (user_id,))
                return (0, "user removed")
            else:
                return (-1, "ERROR: no user with that id!")
        elif 'username' in pardic.keys():
            user_id = self.execute_sql("SELECT id FROM users WHERE username='%s'" % (pardic['username'],))
            if user_id:
                self.execute_sql("DELETE FROM usergroups WHERE user_id='%s'" % (user_id[0][0],))
                self.execute_sql("DELETE FROM users WHERE id='%s';" % (user_id[0][0],))
                return (0, "user removed")
            else:
                return (-1, "ERROR: no user with that email!")
        # TODO: also remove user from a group
        elif 'name' in pardic.keys() or 'group_id' in pardic.keys():
            return (-1, "ERROR: email or id required!")
        else:
            return (-1, "ERROR: User does not exist!")

    def add_group(self, pardic):
        """
        Adds a new group. Checks for duplicates in name. If group exists:

        (-1, "ERROR: Group exists!")

        """
        groups = [des[0] for des in self.execute_sql('SELECT designator FROM groups;')]
        if not pardic['designator']:
            return (-1, 'ERROR: no group designator provided!')
        if pardic['designator'] not in groups:
            self.execute_sql(("INSERT INTO groups (designator) VALUES ('%s')"
                              % (pardic['designator'])))
            return (0, "Group added")
        else:
            return (-1, "ERROR: Group exists!")

    def add_to_group(self, user, group):
        """
        Adds the user as member of the group.
        Checks for duplicates in name. If group exists:

        (-1, "ERROR: Group exists!")

        """
        if user not in [user_id[0] for user_id in self.execute_sql('SELECT id FROM users')]:
            return (-1, "ERROR: user does not exist")
        if group not in [group_id[0] for group_id in self.execute_sql('SELECT id FROM groups')]:
            return (-1, "ERROR: group does not exist")
        usergroups = self.execute_sql('SELECT (user_id, group_id) FROM usergroups')
        if (user, group) in usergroups:
            return (-1, "ERROR: user already in group")
        elif user in [ids[0] for ids in usergroups]:
            return (-1, "ERROR: user in different group")
            # TODO: make it so a user can have multiple groups?
        else:
            self.execute_sql("INSERT INTO usergroups (user_id, group_id) VALUES ('%s','%s')"
                             % (user, group))
            return (0, "User added to group")

    def add_object(self, pardic, objparams):
        """
          Creates a new object in the db with the characterisstics given by the dictionary of parameters: pardic.
          It shall create a new object in periodic or any of the solar system objects (SSO) if necessary.
            The parameters will be specified in the objparams dictionary.
          In case of a SSO, the user has the option to specify only the name if the object is know.
          The function shall check the parameters of the object from the .edb XEphem files and fill the table
            corresponding to the orbital parameters of the object.
          The return of the function can be as follows:

          (CODE, MESSAGE)

          i.e.
          (0, "Object added")
          (-1, "ERROR: the orbital parameters for the SSO are now tabulated. Please, introduce them manually.")
          """
        # TODO: test for each different case of typedesig
        # TODO: have it update existing object if it already exists?
        # TODO: set objparams default to None?

        if pardic['marshal_id'] in [obj[0] for obj in self.execute_sql('SELECT marshal_id FROM object')]:
            return (-1, "ERROR: object exists")
        # TODO: when q3c added to database, search of object withing 1 arcsecond radius (consider anythin inside it as duplicate)
        obj_sql = ("INSERT INTO object (marshal_id, name, ra, dec, typedesig, epoch) Values "
                   "('%s','%s','%s','%s','%s','%s');" % (
                        pardic['marshal_id'], pardic['name'], pardic['ra'], pardic['dec'],
                        pardic['typedesig'], pardic['epoch']))

        self.execute_sql(obj_sql)  # TODO, modify it depending on what parameters make sense for the type
        object_id = self.execute_sql("SELECT id FROM object WHERE marshal_id = %s" % (pardic['marshal_id'],))[0]
        if pardic['typedesig'] == 'f':
            return (0, "Object added")

        elif pardic['typedesig'] == 'e':
            sql = ("INSERT INTO elliptical_heliocentric (object_id, inclination, longascnode_O, perihelion_o, "
                   "a, n, e, M, mjdepoch, D, M1, M2, s) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s',"
                   "'%s','%s','%s','%s')" % (
                    pardic['id'], objparams['inclination'], objparams['longascnode_O'],
                    objparams['perihelion_o'], objparams['a'], objparams['n'], objparams['e'], objparams['M'],
                    objparams['mjdepoch'], objparams['D'], objparams['M1'], objparams['M2'], objparams['s']))
            self.execute_sql(sql)
            return (0, "Elliptical heliocentric object added")

        elif pardic['typedesig'] == 'h':
            sql = (
             "INSERT INTO hyperbolic_heliocentric (object_id, T, inclination, longascnode_O, perihelion_o, "
             "e, q, D, M1, M2, s) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" %
             (pardic['id'], objparams['T'], objparams['inclination'], objparams['longascnode_O'],
              objparams['perihelion_o'], objparams['e'], objparams['q'], objparams['D'], objparams['M1'],
              objparams['M2'], objparams['s']))
            self.execute_sql(sql)
            return (0, "Hyperbolic heliocentric object added")

        elif pardic['typedesig'] == 'p':
            sql = (
             "INSERT INTO parabolic_heliocentric (object_id, T, e, inclination, longascnode_O, perihelion_o, "
             "q, D, M1, M2, s) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" %
             (pardic['id'], objparams['T'], objparams['e'], objparams['inclination'],
              objparams['longascnode_O'], objparams['perihelion_o'], objparams['q'], objparams['D'],
              objparams['M1'], objparams['M2'], objparams['s']))
            self.execute_sql(sql)
            return (0, "Parabolic heliocentric object added")

        elif pardic['typedesig'] == 'E':
            sql = (
             "INSERT INTO earth_satellite (object_id, T, inclination, ra, e, pedigree, M, n, "
             "decay, reforbit, drag) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" %
             (pardic['id'], objparams['T'], objparams['inclination'], objparams['ra'],
              objparams['e'], objparams['pedigree'], objparams['M'], objparams['n'], objparams['decay'],
              objparams['reforbit'], objparams['drag']))
            self.execute_sql(sql)
            return (0, "Earth satellite object added")
            

        elif pardic['typedesig'] == 'P':
            # TODO, use the planet/satellite name (.edb XEphem) to generate the orbit
            obj_name = pardic['name']
            
            pass

    def add_request(self, pardic):
        """
         Creates a new object in the request table with the parameters form the dictionary.
         It shall check that there are no duplicate requests for the same object, obsreving time and project ID.
         If there is a duplicate, it should return a negative result and a message.

         (-1, "ERROR: this request is a duplicate!")
         """
        # TODO: indicate that exptime should be of the form '{ifu_time, rc_time}' with those being time per exposure
        # TODO: remove options that are default null (let them be introduced by update_request)
        # TODO: have update_request handle status, cadence, phasesamples, sampletolerance, odering (,nexposures, filter?)
        requests = self.execute_sql("SELECT (object_id, program_id, marshal_id, exptime, maxairmass,"
               " priority) FROM request WHERE status != 'EXPIRED';")
        # TODO: check mainly for program_id, issue warning that it is repeat, but allow
        if (pardic['object_id'], pardic['program_id'], pardic['marshal_id'],
                pardic['exptime'], pardic['maxairmass'], pardic['priority']) in requests: 
            # warnings.warn
            print "program has already requested this object"
        if pardic['object_id'] not in [obj[0] for obj in self.execute_sql('SELECT id FROM object;')]:
            return (-1, "ERROR: object does not exist")

        sql = ("INSERT INTO request (object_id, user_id, program_id, marshal_id, exptime, maxairmass,"
               " priority, inidate, enddate) VALUES ('%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s');" %
               (pardic['object_id'], pardic['user_id'], pardic['program_id'], pardic['marshal_id'],
                pardic['exptime'], pardic['maxairmass'], pardic['priority'], pardic['inidate'],
                pardic['enddate']))
        self.execute_sql(sql)
        return (0, "Request added")

    def update_request(self, pardic):
        """
         Updates the request table with the parameters from the dictionary.
      
          possible values for status are:
            - PENDING
            - ACTIVE
            - COMPLETED
            - CANCELED
            - EXPIRED

         """
        # TODO: test, complete docstring, restrict which parameters are allowed to be changed?
        keys = pardic.keys()
        keys.remove('id')
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update")
        sql = "UPDATE request SET "
        for key in keys:
            sql += "%s = '%s', " % (key, pardic[key])
        sql += "lastmodified = 'NOW()' "  # TODO: check if this works
        sql += "WHERE id = %s;" % (pardic['id'],)
        # TODO: check if the requests exists yet
        self.execute_sql(sql)
        return (0, "Requests updated")

    def get_active_requests(self):
        """
         Returns the list of requests that are active now.

         """
        # TODO: test, add formatting?
        active_requests = self.execute_sql("SELECT * FROM request WHERE status='ACTIVE';")
        return active_requests

    def expire_requests(self):
        """
         Updates the request table. For all the active requests that were not completed,
          and had an expiry date before than NOW(), are marked as "EXPIRED".

        """
        sql = "UPDATE request SET status = 'EXPIRED' WHERE enddate < 'NOW()' AND status != 'COMPLETED';"
        self.execute_sql(sql)
        return (0, "Requests expired")

    def cancel_scheduled_request(self, requestid):
        """
         Changes the status of the scheduled request to "CANCELED"
         """
        if request_id not in [x[0] for x in self.execute_sql("SELECT id FROM request")]:
            return (-1, "ERROR: request does not exist!")
        # cancel the associated atomicrequests       
        self.execute_sql("UPDATE atomicrequest SET status='CANCELED' WHERE request_id='%s'" % (requestid,))
        self.update_request({'id': requestid, 'status': 'CANCELED'})
        return (0, "Request canceled")

    def update_scheduled_request(self, requestid):
        # TODO: find what this was supposed to be
        """
         Updates the scheduled request with new parameters "CANCELED"
         """

        pass

    def cancel_scheduled_between(self, initime, endtime):
        # TODO: find if this means scheduled or requested
        """
         Cancels all the scheduled requests that were scheduled between the ini and end time.
          Their status should be changed to "CANCELED".
         """

        pass    

    def add_atomic_request(self, pardic):
        """
         Creates a new object in the schedule table with the parameters specified in the dictionary.

        - PENDING
        - OBSERVED
        - REDUCED

         """
        if pardic['request_id'] not in [obj[0] for obj in self.execute_sql('SELECT id FROM request;')]:
            return (-1, "ERROR: request does not exist")
        # TODO: query for object_id from given request_id
        if not pardic['object_id'] == self.execute_sql("SELECT object_id FROM request WHERE id='%s'" % (pardic['request_id']))[0][0]:
            return (-1, "ERROR: object different from request object!")

        sql = ("INSERT INTO request (request_id, user_id, program_id, marshal_id, exptime, maxairmass,"
               " priority, inidate, enddate) VALUES ('%s', '%s', '%s', '%s','%s', '%s', '%s', '%s', '%s');" %
               (pardic['request_id'], pardic['user_id'], pardic['program_id'], pardic['marshal_id'],
                pardic['exptime'], pardic['maxairmass'], pardic['priority'], pardic['inidate'],
                pardic['enddate']))
        self.execute_sql(sql)
        return (0, "Request added")
        pass

    def update_atomic_request(self, pardic):
        """
        Updates an atomic request with the parameters from the dictionary

        """
        # TODO: test, complete docstring, restrict which parameters are allowed to be changed?
        keys = pardic.keys()
        keys.remove('id')
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update")
        sql = "UPDATE atomicrequest SET "
        for key in keys:
            sql += "%s = '%s', " % (key, pardic[key])
        sql += "lastmodified = 'NOW()' "  # TODO: check if this works
        sql += "WHERE id = %s;" % (pardic['id'],)
        self.execute_sql(sql)

    def clear_atomic_request(self, pardic):
        """
        clear the database of atomic requsts at the end of a night
        """
        pass

    def get_request_atomic_requests(self, request_id):
        """
        return the atomicreqests associated with a single request

        Args:
            request_id: int
                The request_id of the desired request

        Returns:
            all atomicreqests associated with the desired request
        """
        # TODO: test
        atomic_requests = self.execute_sql("SELECT * from atomicrequest WHERE request_id = %s" % (request_id,))
        return atomic_requests

    def add_observation_fits(self, fitsfile):
        """
         Receives the fits file path.
         The code shall be able to read the header of the fits file and create an instance in the observation table.

         Whenever this observation is added, the code updates the table "schedule" and marks the scheduled request
         as "observed".


         The code also checks if all the scheduled requests associated with a general request has been observed.
         If it is the case, it shall change the status of the general request to "COMPLETED".

         Finally, using the fits keywords related to telescope status, the function also creates a new register
         in the telescope_stats table.
         """
        hdulist = fits.open(fitsfile)
        header = hdulist[0].header
        #atomic_id = header['ATOM_ID']  # TODO: add ATOM_ID to the header
        atomic_id = 1
        request_id = 1
        object_id =1         
        #request_id, object_id = self.execute_sql("SELECT (request_id, object_id) FROM atomicrequest WHERE id='%s'" % (atomic_id,))[0]
        header_dict = {'object_id': object_id, 'request_id': request_id, 'atomicrequest_id': atomic_id,
                       'mjd': header['JD'], 'airmass': header['AIRMASS'], 'exptime': header['EXPTIME'],
                       'fitsfile': fitsfile, 'lst': header['LST'], 'ra': header['RA'], 'dec': header['DEC'],
                       'tel_ra': header['TEL_RA'], 'tel_dec': header['TEL_DEC'], 'tel_az': header['TEL_AZ'],
                       'tel_el': header['TEL_EL'], 'tel_pa': header['TEL_PA'], 'ra_off': header['RA_OFF'],
                       'dec_off': header['DEC_OFF'], 'utc': header['UTC'], 'camera': header['CAM_NAME']}

        tel_stats = {'date': header['OBSDATE'], 'dome_status': header['DOMEST'], 'in_temp': header['IN_AIR'],
                     'in_humidity': header['IN_HUM'], 'in_dew': header['IN_DEW'], 'out_temp': header['OUT_AIR'],
                     'out_humidity': header['OUT_HUM'], 'out_dew': header['OUT_HUM'], 'wind_dir': header['WIND_DIR'], 
                     'wsp_cur': header['WSP_CUR'], 'wsp_avg': header['WSP_AVG'], 'mir_temp': header['MIR_TEMP'],
                     'top_air': header['TOP_AIR'], 'pri_temp': header['PRI_TEMP'], 'sec_temp': header['SEC_TEMP'],
                     'flo_temp': header['FLO_TEMP'], 'bot_temp': header['BOT_TEMP'], 'mid_temp': header['MID_TEMP'],
                     'top_temp': header['TOP_TEMP']}
        # TODO: add tel_stats['observation_id']
        print header_dict
        print tel_stats
        # TODO: make sure the atomicrequest_id is stored in header
        request_status = [status[0] for status in self.execute_sql("SELECT status FROM atomicrequest WHERE request_id='%s'" % (request_id,))]
        print request_status
        if all(np.array(request_status) == 'OBSERVED'):
            self.update_request({'id': request_id, 'status': 'COMPLETED'})
        pass


    def add_reduced_photometry(self, pardic):
        """
         Creates a new object in the phot table with the parameters specified in the dictionary.
         Only one reduction shall exist for each observation. If the reduction exists, an update
         is made.

          This function also updates the schdele table and changes the status 
          of the associated request to "REDUCED".
         """
        
        pass

    def add_reduced_spectrum(self, pardic):
        """
         Creates a new object in the spec table with the parameters specified in the dictionary.
         Only one reduction shall exist for each observation. If the reduction exists, an update
         is made.

          This function also updates the schdele table and changes the status 
          of the associated request to "REDUCED".
         """
        pass

    def add_metrics_phot(self, pardic):
        """
         Creates a new object in the metrics phot table with the parameters specified in the dictionary.
         Only one metric shall exist for each observation. If the reduction exists, an update
         is made.
         """
        # create and/or update
        pass

    def add_metrics_spec(self, pardic):
        """
         Creates a new object in the metrics spec stats table with the parameters specified in the dictionary.
         Only one metric shall exist for each observation. If the reduction exists, an update
         is made.
         """
        # create table
        pass

    def add_flexure(self, pardic):
        """
         Creates a new object in the flexure table.
         """
        # if flexure is too high, tell something is wrong?
        pass

    def add_classification(self, pardic):
        """
        Creates a classification object attached to the reduced spectrum.
        """
        pass

