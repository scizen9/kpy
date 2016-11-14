import sqlalchemy.pool as pool
from sqlalchemy import exc, create_engine
import psycopg2
import numpy as np
import subprocess
import warnings


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
        # TODO: is name, id, or e-mail the best designator to test for?
        # check if already contained
        users = self.execute_sql('SELECT * FROM users')
        usernames = [user[1] for user in users]
        user_ids = [user[0] for user in users]
        if not (pardic['name'] in usernames or pardic['id'] in user_ids):
            self.execute_sql("INSERT INTO users (id, name, email) Values ('%s', '%s', '%s');"
                             % (pardic['id'], pardic['name'], pardic['email']))
            return (0, "User added")
        elif pardic['name'] in usernames and pardic['id'] in user_ids:
            return (-1, "ERROR: User with that name and id exists")
        elif pardic['name'] in usernames:
            return (-1, "ERROR: User with that name exists!")
        else:
            return (-1, "ERROR: User with that id exists!")
        # TODO: test this

    def remove_user(self, pardic):
        """
        Removes the user. If user does not exist:
          (-1, "ERROR: User does not exist!")

        """
        users = self.execute_sql('SELECT * FROM users')  # [(id, name, email), (id, name, email), ...]
        # TODO: only look at the id?
        if (int(pardic['id']), pardic['name'], pardic['email']) in users:
            self.execute_sql("DELETE FROM users WHERE id='%s';" % (pardic['id'],))
            return (0, "user removed")
        else:
            return (-1, "ERROR: User does not exist!")

    def add_group(self, pardic):
        """
        Adds a new group. Checks for duplicates in name. If group exists:

        (-1, "ERROR: Group exists!")

        """
        # Need a group table or a group_id (in the user table)
        pass

    def add_to_group(self, user, group):
        """
        Adds the user as member of the group.
        Checks for duplicates in name. If group exists:

        (-1, "ERROR: Group exists!")

        """
        pass

    def add_object(self, pardic, objparams):
        """
          Creates a new object in the db with the characterisstics given by the dictionary of parameters: pardic.
          It shall create a new object in periodic or any of the solar system objects (SS) if necessary.
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
        if objparams:
            if not pardic['id'] == objparams['object_id']:
                raise AttributeError('object_id in objparams differs from id in pardic')

        obj_sql = ("INSERT INTO object (id, marshal_id, name, ra, dec, typedesig, epoch) Values "
                   "('%s','%s','%s','%s','%s','%s','%s');" % (
                        pardic['id'], pardic['marshal_id'], pardic['name'], pardic['ra'], pardic['dec'],
                        pardic['typedesig'], pardic['epoch']))

        if not pardic['typedesig'] == 'P':
            self.execute_sql(obj_sql)  # TODO, modify it depending on what parameters make sense for the type
            if pardic['typedesig'] == 'e':
                sql = ("INSERT INTO elliptical_heliocentric (id, object_id, inclination, longascnode_O, perihelion_o, "
                       "a, n, e, M, mjdepoch, D, M1, M2, s) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s',"
                       "'%s','%s','%s','%s')" % (
                        objparams['id'], pardic['id'], objparams['inclination'], objparams['longascnode_O'],
                        objparams['perihelion_o'], objparams['a'], objparams['n'], objparams['e'], objparams['M'],
                        objparams['mjdepoch'], objparams['D'], objparams['M1'], objparams['M2'], objparams['s']))
                self.execute_sql(sql)

            elif pardic['typedesig'] == 'h':
                sql = (
                 "INSERT INTO hyperbolic_heliocentric (id, object_id, T, inclination, longascnode_O, perihelion_o, "
                 "e, q, D, M1, M2, s) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" %
                 (objparams['id'], pardic['id'], objparams['T'], objparams['inclination'], objparams['longascnode_O'],
                  objparams['perihelion_o'], objparams['e'], objparams['q'], objparams['D'], objparams['M1'],
                  objparams['M2'], objparams['s']))
                self.execute_sql(sql)

            elif pardic['typedesig'] == 'p':
                sql = (
                 "INSERT INTO parabolic_heliocentric (id, object_id, T, e, inclination, longascnode_O, perihelion_o, "
                 "q, D, M1, M2, s) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" %
                 (objparams['id'], pardic['id'], objparams['T'], objparams['e'], objparams['inclination'],
                  objparams['longascnode_O'], objparams['perihelion_o'], objparams['q'], objparams['D'],
                  objparams['M1'], objparams['M2'], objparams['s']))
                self.execute_sql(sql)

            elif pardic['typedesig'] == 'E':
                sql = (
                 "INSERT INTO parabolic_heliocentric (id, object_id, T, inclination, ra, e, pedigree, M, n, "
                 "decay, reforbit, drag) VALUES ('%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s','%s')" %
                 (objparams['id'], pardic['id'], objparams['T'], objparams['inclination'], objparams['ra'],
                  objparams['e'], objparams['pedigree'], objparams['M'], objparams['n'], objparams['decay'],
                  objparams['reforbit'], objparams['drag']))
                self.execute_sql(sql)
            else:
                pass

        elif pardic['typedesig'] == 'P':
            # TODO, use the planet/satellite name (.edb XEphem) to generate the orbit
            pass

    def add_request(self, pardic):
        """
         Creates a new object in the request table with the parameters form the dictionary.
         It shall check that there are no duplicate requests for the same object, obsreving time and project ID.
         If there is a duplicate, it should return a negative result and a message.

         (-1, "ERROR: this request is a duplicate!")
         """
        # TODO: remove options that are default null (let them be introduced by update_request)
        # TODO: have update_request handle status, cadence, phasesamples, sampletolerance, odering (,nexposures, filter?)
        sql = ("INSERT INTO request (id, object_id, user_id, program_id, marshal_id, exptime, maxairmass,"
               " priority, inidate, enddate, filters, ellinexposures) VALUES ('%s', '%s', '%s', '%s', '%s',"
               "'%s', '%s', '%s', '%s', '%s', '%s', '%s');" %
               (pardic['id'], pardic['object_id'], pardic['user_id'], pardic['program_id'], pardic['marshal_id'],
                pardic['exptime'], pardic['maxairmass'], pardic['priority'], pardic['inidate'],
                pardic['enddate'], pardic['filters'], pardic['nexposures']))
        self.execute_sql(sql)
        objects = [obj[0] for obj in self.execute_sql("SELECT id FROM object;")]
        if pardic['object_id'] not in objects:  # does the database check this?
            warnings.warn("object associated with request is not in database", UserWarning)
        pass

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

    def expire_requests(self):
        """
         Updates the request table. For all the active requests that were not completed,
          and had an expiry date before than NOW(), are marked as "EXPIRED".

        """

        pass

    def get_active_requests(self):
        """
         Returns the list of requests that are active now.

         """
        # TODO: test, add formatting?
        active_requests = self.execute_sql("SELECT * FROM request WHERE status='ACTIVE';")
        return active_requests

    def add_atomic_request(self, pardic):
        """
         Creates a new object in the schedule table with the parameters specified in the dictionary.

        - PENDING
        - OBSERVED
        - REDUCED

         """

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
         If it is the case, it shall change the status of the general request to "OBSERVED".

         Finally, using the fits keywords related to telescope status, the function also creates a new register
         in the telescope_stats table.
         """

        pass

    def cancel_scheduled_request(self, requestid):
        """
         Changes the status of the scheduled request to "CANCELED"
         """

        pass

    def update_scheduled_request(self, requestid):
        """
         Updates the scheduled request with new parameters "CANCELED"
         """

        pass

    def cancel_scheduled_between(self, initime, endtime):
        """
         Cancels all the scheduled requests that were scheduled between the ini and end time.
          Their status should be changed to "CANCELED".
         """

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
        pass

    def add_metrics_spec(self, pardic):
        """
         Creates a new object in the metrics spec stats table with the parameters specified in the dictionary.
         Only one metric shall exist for each observation. If the reduction exists, an update
         is made.
         """
        pass

    def add_flexure(self, pardic):
        """
         Creates a new object in the flexure table.
         """
        pass

    def add_classification(self, pardic):
        """
        Creates a classification object attached to the reduced spectrum.
        """
        pass

