import sqlalchemy.pool as pool
from sqlalchemy import exc, create_engine
import psycopg2
import numpy as np
import subprocess


# Singleton/SingletonPattern.py

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

    def execute_sql(self, sql):#, ret = True):
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
        """       
        if ret:
            obj = cursor.fetchall()
            conn.close()
            return obj
        else:
            conn.close()
            return None"""
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
        usernames = self.execute_sql("SELECT name FROM users")

        if not pardic['name'] in usernames:
            sql = ("INSERT INTO users (id, name, email) VALUES ('%s', '%s', '%s');"
                             % (pardic['id'],pardic['name'],pardic['email']))
            print sql            
            self.execute_sql(sql)#, False)            
            return (0, "User added")
        else:
            return (-1, "ERROR: User exists!")
        # TODO: test this

    def remove_user(self, pardic):
        """
        Removes the user. If user does not exist:
          (-1, "ERROR: User does not exist!")

        """
        pass

    def add_group(self, pardic):
        """
        Adds a new group. Checks for duplicates in name. If group exists:

        (-1, "ERROR: Group exists!")

        """
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
          The function shall check the parameters of the object from the .edb XEphem files and fill the table corresponding
            to the orbital parameters of the object.
          The return of the function can be as follows:

          (CODE, MESSAGE)

          i.e.
          (0, "Object added")
          (-1, "ERROR: the orbital parameters for the SSO are now tabulated. Please, introduce them manually.")
          """

        pass

    def add_request(self, pardic):
        """
         Creates a new object in the request table with the parameters form the dictionary.
         It shall check that there are no duplicate requests for the same object, obsreving time and project ID.
         If there is a duplicate, it should return a negative result and a message.

         (-1, "ERROR: this request is a duplicate!")
         """

        pass

    def update_request(self, pardic):
        """
         Updates the request table with the parameters form the dictionary.
      
          possible values for status are:
            - PENDIG
            - ACTIVE
            - COMPLETED
            - CANCELED
            - EXPIRED

         """

        pass

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

        pass

    def add_atomic_request(self, pardic):
        """
         Creates a new object in the schedule table with the parameters specified in the dictionary.

        - PENDING
        - OBSERVED
        - REDUCED

         """

        pass

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
