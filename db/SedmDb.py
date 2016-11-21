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
        adds a new user
        Args:
            pardic: dict
                required: 'username' (unique), 'name', 'email'
        Returns:
            (-1, "ERROR: Username exists") if the username is a duplicate
            (0, "....") if the user was added
        """

        # check if already contained
        usernames = [user[0] for user in self.execute_sql('SELECT username FROM users')]
        if not pardic['username'] in usernames:
            self.execute_sql("INSERT INTO users (username, name, email) Values ('%s', '%s', '%s');"
                             % (pardic['username'], pardic['name'], pardic['email']))
            user_id = self.execute_sql("SELECT id FROM users WHERE username='%s'" % (pardic['username'],))[0][0]

            # TODO: choose between group_id and designator or handle both?
            if 'group_id' in pardic.keys():
                group_add = self.add_to_group(user_id, pardic['group_id'])
            elif 'group_designator' in pardic.keys():
                group_id = self.execute_sql("SELECT id FROM groups WHERE designator='%s'"
                                            % (pardic['group_designator'],))[0][0]
                group_add = self.add_to_group(user_id, group_id)
            else:
                group_add = (1, "no group provided")

            if group_add[0] == -1:
                return (0, "User added, failed to add to group")
            elif group_add[0] == 1:
                return (0, "User added, no group provided")
            else:
                return (0, "User added")
        else:
            return (-1, "ERROR: User with that username exists!")
        # TODO: test

    def remove_user(self, pardic):
        """

        Args:
            pardic: dict
                required: 'username' (recommended) OR 'id'

        Returns:
            (-1, "ERROR: ...") if there was an issue (user doesn't exist, not enough information in pardic)
            (0, "User removed") if the removal was successful
        """
        if 'username' in pardic.keys():
            user_id = self.execute_sql("SELECT id FROM users WHERE username='%s'" % (pardic['username'],))
            if user_id:
                self.execute_sql("DELETE FROM usergroups WHERE user_id='%s'" % (user_id[0][0],))
                self.execute_sql("DELETE FROM users WHERE id='%s';" % (user_id[0][0],))
                return (0, "User removed")
            else:
                return (-1, "ERROR: no user with that username!")
        elif 'id' in pardic.keys():
            if pardic['id'] in [x[0] for x in self.execute_sql('SELECT id FROM users;')]:
                self.execute_sql("DELETE FROM usergroups WHERE id='%s'" % (pardic['id'],))
                self.execute_sql("DELETE FROM users WHERE id='%s';" % (pardic['id'],))
                return (0, "User removed")
            else:
                return (-1, "ERROR: no user with that id!")
        elif 'name' in pardic.keys() or 'email' in pardic.keys():
            return (-1, "ERROR: username or id required!")
        else:
            return (-1, "ERROR: User does not exist!")
        # TODO: test id, username, neither and incorrect ones

    def add_group(self, pardic):
        """
        Adds a new group. Checks for duplicates in name
        Args:
            pardic: dict
                required: 'designator'

        Returns:
            (-1, "ERROR: ...") if no designator was provided or there is already a group with it
            (0, "Group added") if the adding was successful
        """
        if not pardic['designator']:
            return (-1, 'ERROR: no group designator provided!')
        groups = [des[0] for des in self.execute_sql('SELECT designator FROM groups;')]
        if pardic['designator'] not in groups:
            self.execute_sql(("INSERT INTO groups (designator) VALUES ('%s')"
                              % (pardic['designator'])))
            return (0, "Group added")
        else:
            return (-1, "ERROR: Group exists!")

    def add_to_group(self, user, group):
        """
        Adds the user as member of the group. Checks for duplicates in name.
        Args:
            user: int
                id of the user in the 'users' Table
            group: int
                id of the group in the 'groups' Table

        Returns:
            (-1, "ERROR: ...") if there was areason for failure
            (0, "User added to group") if the adding was successful
        """
        if user not in [user_id[0] for user_id in self.execute_sql('SELECT id FROM users')]:
            return (-1, "ERROR: user does not exist")
        if group not in [group_id[0] for group_id in self.execute_sql('SELECT id FROM groups')]:
            return (-1, "ERROR: group does not exist")
        usergroups = self.execute_sql('SELECT (user_id, group_id) FROM usergroups')
        if (user, group) in usergroups:
            return (-1, "ERROR: user already in group")
        else:
            self.execute_sql("INSERT INTO usergroups (user_id, group_id) VALUES ('%s','%s')"
                             % (user, group))
            return (0, "User added to group")
        # TODO: test adding a user to group, non-existant user, non-existant group, already connected user/group

    def add_object(self, pardic, objparams):
        """
        Creates a new object, if the object is a solar system object (SSO) it attempts to find an entry for it
        in the .edb (XEphem) files. If objparams are given, they are used (and update the file if the object is found)
        and are stored in the corresponding Table for the objects orbit type.
        Args:
            pardic: dict
                requred: 'name', 'typedesig'(, 'ra', 'dec', 'epoch' OR 'marshal_id' for fixed object)
                optional: 'iauname', 'marshal_id', 'epoch', 'ra', 'dec'
                NOTE: 'typedesig' should be one of:
                            'f' (fixed), 'P' (built-in planet or satellite name), 'e' (heliocentric elliptical),
                            'h' (heliocentric hyperbolic), 'p' (heliocentric parabolic), 'E' (geocentric elliptical)
            objparams: dict (different required parameters needed for each typedesig)
                'f' or 'P' (or if pardic['name'] is in .edb file): none needed
                'e': 'inclination', 'longascnode_0' (lon. of ascending node), 'perihelion_o' (arg. of perihelion),
                     'a' (mean distance AU), 'n' (mean daily motion deg/day), 'e' (eccentricity),
                     'M' (mean anomaly), 'mjdepoch' (epoch, time of 'M'), 'D' (equinox year),
                     'M1', 'M2' (first and second components of magnitude model), 's' (angular size at 1 AU)
                'h': 'T' (date), 'inclination', 'longascnode_0' (lon. of ascending node),
                     'perihelion_o' (arg. of perihelion), 'e' (eccentricity), 'q' (perihelion distance AU),
                     'D' (equinox year), 'M1', 'M2' (first and second components of magnitude model),
                     's' (angular size at 1 AU)
                'p': 'T' (date), 'inclination', 'longascnode_0' (lon. of ascending node),
                     'perihelion_o' (arg. of perihelion), 'q' (perihelion distance), 'D' (equinox year),
                     'M1', 'M2' (first and second components of magnitude model), 's' (angular size at 1 AU)
                'E': 'T' (date), 'inclination', 'ra' (ra of ascending node), 'e' (eccentricity),
                     'pedigree' (arg. of pedigree), 'M' (mean anomaly), 'n' (mean motion, revs/day),
                     'decay' (orbit decay rate, rev/day^2), 'reforbit' (integral reference orbit number at epoch),
                     'drag' (drag coefficient, 1/(Earth radii))


        Returns:
            (-1, "ERROR: ...") if it failed to add
            (0, "Object added") if the object is added successfully
        """


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
        # TODO: for each typedesig, handle "name" and search the .edb for it (instead of requiring objparams)
        # TODO: test for each different case of typedesig
        # TODO: have it update existing object if it already exists?
        # TODO: set objparams default to None?
        # TODO: add iauname to all of the obj_sql and generate them rather than hard-code
        if 'marshal_id' in pardic.keys():
            if pardic['marshal_id'] in [obj[0] for obj in self.execute_sql('SELECT marshal_id FROM object')]:
                return (-1, "ERROR: object exists")
        # TODO: when q3c added to database, search of object withing 1 arcsecond radius (consider anything inside it as duplicate)
        if pardic['typedesig'] == 'f':
            obj_sql = ("INSERT INTO object (marshal_id, name, ra, dec, typedesig, epoch) Values "
                       "('%s','%s','%s','%s','%s','%s');" % (
                           pardic['marshal_id'], pardic['name'], pardic['ra'], pardic['dec'],
                           pardic['typedesig'], pardic['epoch']))

            self.execute_sql(obj_sql)
            return (0, "Object added")

        elif pardic['typedesig'] == 'e':
            # TODO: modify obj_sql for the specific typedesig
            obj_sql = ("INSERT INTO object (marshal_id, name, ra, dec, typedesig, epoch) Values "
                       "('%s','%s','%s','%s','%s','%s');" % (
                           pardic['marshal_id'], pardic['name'], pardic['ra'], pardic['dec'],
                           pardic['typedesig'], pardic['epoch']))
            self.execute_sql(obj_sql)
            object_id = self.execute_sql("SELECT id FROM object WHERE marshal_id = %s" % (pardic['marshal_id'],))[0]

            sql = ("INSERT INTO elliptical_heliocentric (object_id, inclination, longascnode_O, perihelion_o, "
                   "a, n, e, M, mjdepoch, D, M1, M2, s) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s',"
                   "'%s','%s','%s','%s')" % (
                    object_id, objparams['inclination'], objparams['longascnode_O'],
                    objparams['perihelion_o'], objparams['a'], objparams['n'], objparams['e'], objparams['M'],
                    objparams['mjdepoch'], objparams['D'], objparams['M1'], objparams['M2'], objparams['s']))
            self.execute_sql(sql)
            return (0, "Elliptical heliocentric object added")

        elif pardic['typedesig'] == 'h':
            # TODO: modify obj_sql for the specific typedesig
            obj_sql = ("INSERT INTO object (marshal_id, name, ra, dec, typedesig, epoch) Values "
                       "('%s','%s','%s','%s','%s','%s');" % (
                           pardic['marshal_id'], pardic['name'], pardic['ra'], pardic['dec'],
                           pardic['typedesig'], pardic['epoch']))
            self.execute_sql(obj_sql)
            object_id = self.execute_sql("SELECT id FROM object WHERE marshal_id = %s" % (pardic['marshal_id'],))[0]

            sql = (
             "INSERT INTO hyperbolic_heliocentric (object_id, T, inclination, longascnode_O, perihelion_o, "
             "e, q, D, M1, M2, s) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" %
             (object_id, objparams['T'], objparams['inclination'], objparams['longascnode_O'],
              objparams['perihelion_o'], objparams['e'], objparams['q'], objparams['D'], objparams['M1'],
              objparams['M2'], objparams['s']))
            self.execute_sql(sql)
            return (0, "Hyperbolic heliocentric object added")

        elif pardic['typedesig'] == 'p':
            # TODO: modify obj_sql for the specific typedesig
            obj_sql = ("INSERT INTO object (marshal_id, name, ra, dec, typedesig, epoch) Values "
                       "('%s','%s','%s','%s','%s','%s');" % (
                           pardic['marshal_id'], pardic['name'], pardic['ra'], pardic['dec'],
                           pardic['typedesig'], pardic['epoch']))
            self.execute_sql(obj_sql)
            object_id = self.execute_sql("SELECT id FROM object WHERE marshal_id = %s" % (pardic['marshal_id'],))[0]

            sql = (
             "INSERT INTO parabolic_heliocentric (object_id, T, e, inclination, longascnode_O, perihelion_o, "
             "q, D, M1, M2, s) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" %
             (object_id, objparams['T'], objparams['e'], objparams['inclination'],
              objparams['longascnode_O'], objparams['perihelion_o'], objparams['q'], objparams['D'],
              objparams['M1'], objparams['M2'], objparams['s']))
            self.execute_sql(sql)
            return (0, "Parabolic heliocentric object added")

        elif pardic['typedesig'] == 'E':
            # TODO: modify obj_sql for the specific typedesig
            obj_sql = ("INSERT INTO object (marshal_id, name, ra, dec, typedesig, epoch) Values "
                       "('%s','%s','%s','%s','%s','%s');" % (
                           pardic['marshal_id'], pardic['name'], pardic['ra'], pardic['dec'],
                           pardic['typedesig'], pardic['epoch']))
            self.execute_sql(obj_sql)
            object_id = self.execute_sql("SELECT id FROM object WHERE marshal_id = %s" % (pardic['marshal_id'],))[0]

            sql = (
             "INSERT INTO earth_satellite (object_id, T, inclination, ra, e, pedigree, M, n, "
             "decay, reforbit, drag) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" %
             (object_id, objparams['T'], objparams['inclination'], objparams['ra'],
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

        Args:
            pardic: dict
                required: 'object_id', 'user_id', 'program_id', 'exptime' (string '{spec_duration, phot_duration}'),
                          'priority', 'inidate' (start of observing window), 'enddate' (end of observing window)
                optional: 'marshal_id', 'maxairmass' (max allowable airmass for observation, default 2.5),
                          'cadence' (


        Returns:

        """


        """
         Creates a new object in the request table with the parameters form the dictionary.
         It shall check that there are no duplicate requests for the same object, obsreving time and project ID.
         If there is a duplicate, it should return a negative result and a message.

         (-1, "ERROR: this request is a duplicate!")
         """
        # TODO: handle exptime in-function?
        # TODO: indicate that exptime should be of the form '{ifu_time, rc_time}' with those being time per exposure
        requests = self.execute_sql("SELECT (object_id, program_id) FROM request WHERE status != 'EXPIRED';")
        # TODO: check mainly for program_id, issue warning that it is a repeat, but allow
        if (pardic['object_id'], pardic['program_id']) in requests:
            # warnings.warn
            print "program %s has already requested object: %s" % (pardic['program_id'], pardic['object_id'])
        if pardic['object_id'] not in [obj[0] for obj in self.execute_sql('SELECT id FROM object;')]:
            return (-1, "ERROR: object does not exist")
        # TODO: set default inidate/enddate?

        default_params = ['object_id', 'user_id', 'program_id', 'exptime', 'priority',
                          'inidate', 'enddate']
        columns = '('
        values = "("
        # add the required parameters
        for param in default_params:
            if param not in pardic.keys():
                return (-1, "ERROR: %s not in dictionary" % (param,))
            elif pardic[param]:
                columns += (param + ", ")  # e.g. "(cadence, sampletolerance, "
                values += "'%s', " % (pardic[param],)  # e.g. "('24000', '.2', "
            else:
                return (-1, "ERROR: no value provided for %s" % (param,))
        # add non-required parameters
        other_params = [param for param in pardic.keys() if param not in default_params]
        for param in other_params:
            # if there is actually a values associated with it
            if pardic[param]:
                columns += (param + ", ")
                values += "'%s', " % (pardic[param],)
        columns = columns[:-2] + ')'
        values = values[:-2] + ')'
        sql = ("INSERT INTO request %s VALUES %s;" % (columns, values))
        self.execute_sql(sql)
        return (0, "Request added")
        # TODO: test ...

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
        keys = list(pardic.keys())
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
        # TODO: test...

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
        # tests written
        sql = "UPDATE request SET status = 'EXPIRED' WHERE enddate < 'NOW()' AND status != 'COMPLETED';"
        self.execute_sql(sql)
        return (0, "Requests expired")

    def cancel_scheduled_request(self, requestid):
        """
         Changes the status of the scheduled request to "CANCELED"
         """
        if requestid not in [x[0] for x in self.execute_sql("SELECT id FROM request")]:
            return (-1, "ERROR: request does not exist!")
        # cancel the associated atomicrequests       
        self.execute_sql("UPDATE atomicrequest SET status='CANCELED' WHERE request_id='%s'" % (requestid,))
        self.update_request({'id': requestid, 'status': 'CANCELED'})
        return (0, "Request canceled")
        # TODO: write tests for existing, non-existing request

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

        Args:
            pardic:

                filter options: 'u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b'
        Returns:

        """


        """
         Creates a new object in the schedule table with the parameters specified in the dictionary.

        - PENDING
        - OBSERVED
        - REDUCED

         """
        req_obj_stat = self.execute_sql("SELECT (object_id, status) FROM request WHERE id='%s'" % (pardic['request_id']))
        if not req_obj_stat:  # if there is no request with the id given
            return (-1, "ERROR: request does not exist!")
        # TODO: make sure docstring indicates object_id, ... are generated
        pardic['object_id'] = req_obj_stat[0]  # set object_id given the request
        if req_obj_stat[1] == 'EXPIRED':
            return (-1, "ERROR: request has expired!")

        sql = ("INSERT INTO request (object_id, request_id, order_id, exptime, filter, maxairmass,"
               " priority, inidate, enddate) VALUES ('%s', '%s', '%s','%s', '%s', '%s', '%s', '%s');" %
               (pardic['object_id'], pardic['request_id'], pardic['order_id'],
                pardic['exptime'],  pardic['filter'], pardic['priority'],
                pardic['inidate'], pardic['enddate']))
        # atomic_requests should be created at the start, and purged at the end of each night
        # TODO: given above comment, do we need inidate/enddate?
        self.execute_sql(sql)
        return (0, "Request added")
        # TODO: test

    def update_atomic_request(self, pardic):
        """
        Updates an atomic request with the parameters from the dictionary

        """
        # TODO: test, complete docstring, restrict which parameters are allowed to be changed?
        keys = list(pardic.keys())
        keys.remove('id')
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update")
        sql = "UPDATE atomicrequest SET "
        for key in keys:
            sql += "%s = '%s', " % (key, pardic[key])
        sql += "lastmodified = 'NOW()' "  # TODO: check if this works
        sql += "WHERE id = %s;" % (pardic['id'],)
        self.execute_sql(sql)

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
        # atomic_id = header['ATOM_ID']  # TODO: add ATOM_ID to the header
        atomic_id = 1
        request_id = 1
        object_id = 1
        # request_id, object_id = self.execute_sql("SELECT (request_id, object_id) FROM atomicrequest WHERE id='%s'" % (atomic_id,))[0]
        # TODO: uncomment above and remove the "request_id=1 and object_id=1" once implemented
        header_dict = {'object_id': object_id, 'request_id': request_id, 'atomicrequest_id': atomic_id,
                       'mjd': header['JD'], 'airmass': header['AIRMASS'], 'exptime': header['EXPTIME'],
                       'fitsfile': fitsfile, 'lst': header['LST'], 'ra': header['RA'], 'dec': header['DEC'],
                       'tel_ra': header['TEL_RA'], 'tel_dec': header['TEL_DEC'], 'tel_az': header['TEL_AZ'],
                       'tel_el': header['TEL_EL'], 'tel_pa': header['TEL_PA'], 'ra_off': header['RA_OFF'],
                       'dec_off': header['DEC_OFF'], 'utc': header['UTC'], 'camera': header['CAM_NAME']}
        # generate string describing the columns
        columns = list(header_dict.keys())
        cols = "("
        for col in columns:
            cols += "%s, " % (col,)
        cols = cols[:-2] + ')'  # get rid of extra ", "
        # need to do same with values to have the same parameter order
        values_list = list(header_dict.values())
        values = '(' + str(values_list)[1:-1] + ')'  # make it into the correct value format
        sql = ("INSERT INTO observation %s VALUES %s" % (cols, values))
        self.execute_sql(sql)

        # TODO: see if OUTPUT is able to be used instead of the following
        observation_id = (self.execute_sql("SELECT id FROM observation WHERE atomicrequest_id = '%s'")
                          % (header_dict['atomicrequest_id'],))[0][0]

        tel_stats = {'date': header['OBSDATE'], 'dome_status': header['DOMEST'], 'in_temp': header['IN_AIR'],
                     'in_humidity': header['IN_HUM'], 'in_dew': header['IN_DEW'], 'out_temp': header['OUT_AIR'],
                     'out_humidity': header['OUT_HUM'], 'out_dew': header['OUT_HUM'], 'wind_dir': header['WIND_DIR'], 
                     'wsp_cur': header['WSP_CUR'], 'wsp_avg': header['WSP_AVG'], 'mir_temp': header['MIR_TEMP'],
                     'top_air': header['TOP_AIR'], 'pri_temp': header['PRI_TEMP'], 'sec_temp': header['SEC_TEMP'],
                     'flo_temp': header['FLO_TEMP'], 'bot_temp': header['BOT_TEMP'], 'mid_temp': header['MID_TEMP'],
                     'top_temp': header['TOP_TEMP'], 'observation_id': observation_id}
        # perform the same sql-creation as above
        stat_columns = list(tel_stats.keys())
        stat_cols = "("
        for col in stat_columns:
            stat_cols += "%s, " % (col,)
        stat_cols = stat_cols[:-2] + ')'  # get rid of extra ", "
        stat_values = '(' + str(list(tel_stats.values()))[1:-1] + ')'  # make it into the correct value format
        stat_sql = ("INSERT INTO telescope_stats %s VALUES %s" % (stat_cols, stat_values))
        self.execute_sql(stat_sql)

        # TODO: make sure the atomicrequest_id is stored in header
        request_status = [status[0] for status in
                          self.execute_sql("SELECT status FROM atomicrequest WHERE request_id='%s'" % (request_id,))]
        print request_status  # still testing if this is working
        if all(np.array(request_status) == 'OBSERVED'):
            self.update_request({'id': request_id, 'status': 'COMPLETED'})
        # TODO: add other returns for failure cases?
        return (0, "Observation added")
        # TODO: test with actual fits files

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

