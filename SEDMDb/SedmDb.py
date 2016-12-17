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

        self.sso_objects = None

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
        if 'SELECT' in sql[:8]:
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

        Parameters
        ----------
        pardic: dict
            required: 'username' (unique), 'name', 'email'

        Returns
        -------
        (-1, "ERROR: Username exists") if the username is a duplicate
        (0, "User added") if the user was added
        """
        keys = list(pardic.keys())
        if 'username' not in keys:
            return (-1, "ERROR: no username provided!")
        # check for duplicate username
        usernames = [user[0] for user in self.execute_sql('SELECT username FROM users')]
        if pardic['username'] in usernames:
            return (-1, "ERROR: user with that username exists!")
        for key in reversed(keys):  # remove group keys and any other bad keys
            if key not in ['username', 'name', 'email']:
                keys.remove(key)
        sql = generate_insert_sql(pardic, keys, 'users')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_user sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_user sql command failed with a ProgrammingError!")
        return (0, "User added")

    def remove_user(self, pardic):
        """
        Removes an existing user

        Parameters
        ----------
        pardic: dict
            required: 'username' (recommended) OR 'id'

        Returns
        -------
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
        else:
            return (-1, "ERROR: username or id required!")

    def select_from_users(self, values, where_dict):
        """
        select values from `users`

        Parameters
        ----------
        values: list of str
            values to be returned
        where_dict: dict
            'param':'value' to be used as WHERE clauses
        values/keys options: ['id', 'username', 'name', 'email']

        Returns
        -------
        list: tuples containing the values for each user matching the criteria, empty if no results
        """
        # TODO: test, reconsider return styles
        allowed_params = ['id', 'username', 'name', 'email']
        sql = generate_select_sql(values, where_dict, allowed_params, 'users')

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_group(self, pardic):
        """
        Adds a new group. Checks for duplicates in name

        Parameters
        ----------
        pardic: dict
            required: 'designator'

        Returns
        -------
        (-1, "ERROR: ...") if no designator was provided or there is already a group with it
        (0, "Group added") if the adding was successful
        """
        if 'designator' not in pardic.keys():
            return (-1, 'ERROR: no group designator provided!')
        groups = [des[0] for des in self.execute_sql('SELECT designator FROM groups;')]
        if pardic['designator'] not in groups:
            sql = ("INSERT INTO groups (designator) VALUES ('%s')" % (pardic['designator']))
            try:
                self.execute_sql(sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_group sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_group sql command failed with a ProgrammingError!")
            return (0, "Group added")
        else:
            return (-1, "ERROR: group exists!")

    def add_to_group(self, user, group):
        """
        Adds the user as member of the group. Checks for duplicates in name.

        Parameters
        ----------
        user: int
            id of the user in the 'users' Table
        group: int
            id of the group in the 'groups' Table

        Returns
        -------
        (-1, "ERROR: ...") if there was areason for failure
        (0, "User added to group") if the adding was successful
        """
        if user not in [user_id[0] for user_id in self.execute_sql('SELECT id FROM users')]:
            return (-1, "ERROR: user does not exist!")
        if group not in [group_id[0] for group_id in self.execute_sql('SELECT id FROM groups')]:
            return (-1, "ERROR: group does not exist!")
        usergroups = self.execute_sql('SELECT user_id, group_id FROM usergroups')
        if (user, group) in usergroups:
            return (-1, "ERROR: user already in group!")
        else:
            sql = "INSERT INTO usergroups (user_id, group_id) VALUES ('%s','%s')" % (user, group)
            try:
                self.execute_sql(sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_to_group sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_to_group sql command failed with a ProgrammingError!")
            return (0, "User added to group")

    def add_object(self, pardic):
        """
        Creates a new object

        Parameters
        ----------
        pardic: dict
            required:
                'name', 'typedesig'(, 'ra', 'dec', 'epoch' for a fixed object)
            optional:
                'iauname', 'marshal_id', 'epoch', 'ra', 'dec'
            'typedesig' should be one of:
                'f' (fixed), 'P' (built-in planet or satellite name), 'e' (heliocentric elliptical),
                'h' (heliocentric hyperbolic), 'p' (heliocentric parabolic), 'E' (geocentric elliptical)

        Returns
        -------
        (-1, "ERROR: ...") if it failed to add
        (0, "Object added") if the object is added successfully
        """
        # TODO: have it update existing object if it already exists?
        if 'marshal_id' in pardic.keys():
            if pardic['marshal_id'] in [obj[0] for obj in self.execute_sql('SELECT marshal_id FROM object')]:
                return (-1, "ERROR: object exists!")
        # TODO: when q3c added to database, search of object withing 1 arcsecond radius (consider anything inside it as duplicate)
        obj_keys = list(pardic.keys())

        for key in ['name', 'typedesig']:  # check if 'name' and 'typedesig' are provided
            if key not in obj_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
            elif not pardic[key]:
                return (-1, "ERROR: no value provided for %s!" % (key,))
        for key in reversed(obj_keys):  # remove any extraneous keys
            if key not in ['name', 'typedesig', 'ra', 'dec', 'epoch', 'marshal_id', 'iauname']:
                obj_keys.remove(key)
        pardic['name'] = pardic['name'].lower()  # make all of the names the same format for consistant searching
        if pardic['typedesig'] == 'f':
            for key in ['ra', 'dec', 'epoch']:
                if key not in obj_keys:
                    return (-1, "ERROR: %s not provided!" % (key,))
                elif not pardic[key]:
                    return(-1, "ERROR: no value provided for %s!" % (key,))

            obj_sql = generate_insert_sql(pardic, obj_keys, 'object')
            try:
                self.execute_sql(obj_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_object sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_object sql command failed with a ProgrammingError!")
            return (0, "Fixedd object added")
        else:
            obj_sql = generate_insert_sql(pardic, obj_keys, 'object')
            try:
                self.execute_sql(obj_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_object sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_object sql command failed with a ProgrammingError!")

#        elif not pardic['iauname'] and not orbit_params:  # if not fixed or default, need identifier
#            return (-1, "ERROR: need iauname or orbit_params for non-fixed objects!")
#        elif not orbit_params:
#            return (-1, "ERROR: generating orbit_params from iauname not yet implemented")

    def select_from_objects(self, values, where_dict):
        """
        select values from `objects`

        Parameters
        ----------
        values: list of str
            values to be returned
        where_dict: dict
            'param':'value' to be used as WHERE clauses
        values/keys options: ['id', 'marshal_id', 'name', 'iauname', 'ra', 'dec', 'typedesig', 'epoch']

        Returns
        -------
        list: tuples containing the values for each user matching the criteria, empty if no results
        """
        # TODO: test, reconsider return styles
        allowed_params = ['id', 'marshal_id', 'name', 'iauname', 'ra', 'dec', 'typedesig', 'epoch']
        sql = generate_select_sql(values, where_dict, allowed_params, 'object')

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

# TODO: have get_object for ra/dec
    def get_object_id_from_name(self, object_name):
        """
        finds the id of an object given its name or part of its name

        Parameters
        ----------
        object_name: str

        Returns
        -------
        id, full name if one object is found
        list of all (id, full name) if multiple matches are found for the name
        None if the object is not found
        """
        object_name = object_name.lower()
        sql = "SELECT id, name FROM object WHERE name LIKE '%s%s%s'" % ('%', object_name, '%')
        obj = self.execute_sql(sql)
        if not obj:
            return None
        elif len(obj) > 1:
            return obj
        else:
            return obj[0]  # the sql returns ((id, name),)

# TODO: move following below add_object
# TODO: query associated table for object already existing

    def add_elliptical_orbit(self, orbit_params):
        """
        Adds the orbit parameters for an object

        Parameters
        ----------
        orbit_params: dict
            required:
                'object_id', 'inclination', 'longascnode_0' (lon. of ascending node),
                'perihelion_o' (arg. of perihelion), 'a' (mean distance AU), 'n' (mean daily motion deg/day),
                'e' (eccentricity), 'M' (mean anomaly), 'mjdepoch' (epoch, time of 'M'), 'D' (equinox year),
                'M1', 'M2' (first and second components of magnitude model)

            optional: 's' (angular size at 1 AU)

        Returns
        -------

        """
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['inclination', 'longascnode_0', 'perihelion_o', 'a', 'n', 'e',
                    'M', 'mjdepoch', 'D', 'M1', 'M2', 'object_id']:
            if key not in orb_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
            elif not orbit_params[key]:
                return (-1, "ERROR: no value provided for %s!" % (key,))
        for key in reversed(orb_keys):
            if key not in ['inclination', 'longascnode_0', 'perihelion_o', 'a', 'n', 'e',
                           'M', 'mjdepoch', 'D', 'M1', 'M2', 's', 'object_id']:
                orb_keys.remove(key)

        orb_sql = generate_insert_sql(orbit_params, orb_keys, 'elliptical_heliocentric')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_elliptical_orbit sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_elliptical_orbit sql command failed with a ProgrammingError!")
        return (0, "Elliptical heliocentric orbit added")

    def add_parabolic_orbit(self, orbit_params):
        """


        Parameters
        ----------
        orbit_params: dict
            required:
                'object_id', 'T' (date), 'inclination', 'longascnode_0' (lon. of ascending node),
                'perihelion_o' (arg. of perihelion), 'q' (perihelion distance), 'D' (equinox year),
                'M1', 'M2' (first and second components of magnitude model)

            optional: 's' (angular size at 1 AU)

        Returns
        -------

        """
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['T', 'inclination', 'longascnode_0', 'perihelion_o', 'e', 'q', 'D',
                    'M1', 'M2', 'object_id']:
            if key not in orb_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
            elif not orbit_params[key]:
                return (-1, "ERROR: no value provided for %s!" % (key,))
        for key in reversed(orb_keys):
            if key not in ['T', 'inclination', 'longascnode_0', 'perihelion_o', 'e', 'q', 'D',
                           'M1', 'M2', 's', 'object_id']:
                orb_keys.remove(key)
        orb_sql = generate_insert_sql(orbit_params, orb_keys, 'parabolic_heliocentric')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_parabolic_orbit sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_parabolic_orbit sql command failed with a ProgrammingError!")
        return (0, "Parabolic heliocentric orbit added")

    def add_hyperbolic_orbit(self, orbit_params):
        """

        Parameters
        ----------
        orbit_params: dict
            required:
                'object_id', 'T' (date), 'inclination', 'longascnode_0' (lon. of ascending node),
                'perihelion_o' (arg. of perihelion), 'e' (eccentricity), 'q' (perihelion distance AU),
                'D' (equinox year), 'M1', 'M2' (first and second components of magnitude model)

            optional: 's' (angular size at 1 AU)

        Returns
        -------
        """
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['T', 'inclination', 'longascnode_0', 'perihelion_o', 'e', 'q', 'D',
                    'M1', 'M2', 'object_id']:
            if key not in orb_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
            elif not orbit_params[key]:
                return (-1, "ERROR: no value provided for %s!" % (key,))
        for key in reversed(orb_keys):
            if key not in ['T', 'inclination', 'longascnode_0', 'perihelion_o', 'e', 'q', 'D',
                           'M1', 'M2', 's', 'object_id']:
                orb_keys.remove(key)
        orb_sql = generate_insert_sql(orbit_params, orb_keys, 'hyperbolic_heliocentric')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_hyperbolic_orbit sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_hyperbolic_orbit sql command failed with a ProgrammingError!")
        return (0, "Hyperbolic heliocentric orbit added")

    def add_earth_satellite_orbit(self, orbit_params):
        """

        Parameters
        ----------
        orbit_params: dict
            required:
                'object_id', ' T' (epoch of other fields), 'inclination', 'ra' (ra of ascending node),
                'e' (eccentricity), 'pedigree' (arg. of pedigree), 'M' (mean anomaly),
                'n' (mean motion, revs/day), 'decay' (orbit decay rate, rev/day^2),
                'reforbit' (integral reference orbit number at epoch),

            optional: 'drag' (drag coefficient, 1/(Earth radii))

        Returns
        -------

        """
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['T', 'inclination', 'ra', 'e', 'pedigree', 'M', 'n',
                    'decay', 'reforbit', 'object_id']:
            if key not in orb_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
            elif not orbit_params[key]:
                return (-1, "ERROR: no value provided for %s!" % (key,))
        for key in reversed(orb_keys):
            if key not in ['T', 'inclination', 'ra', 'e', 'pedigree', 'M', 'n',
                           'decay', 'reforbit', 'drag', 'object_id']:
                orb_keys.remove(key)
        orb_sql = generate_insert_sql(orbit_params, orb_keys, 'earth_satellite')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_earth_satellite_orbit sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_earth_satellite_orbit sql command failed with a ProgrammingError!")
        return (0, "Earth satellite orbit added")

    def add_planet_satellite_orbit(self, orbit_params):
        raise NotImplementedError
        # TODO: actually implement? maybe let higher level do this
        # TODO: query ??? for if it already exists

    def add_request(self, pardic):
        """
        Add a request

        Parameters
        ----------
        pardic: dict
            required:
                'object_id', 'user_id', 'program_id', 'exptime' (string '{spec_duration, phot_duration}'),
                'priority', 'inidate' (start of observing window), 'enddate' (end of observing window),
                (one of 'nexposures' or 'ordering' described below)

            optional:
                'marshal_id', 'maxairmass' (max allowable airmass for observation, default 2.5),
                'cadence' (,
                'phasesamples' (how many samples in a period), 'sampletolerance' (how much tolerance in when
                the samples can be taken), 'nexposures' (string '{# of ifu, # of u, # of g, # of r, # of i}'),
                'ordering' (string e.g. '{3g, 3r, 1i, 1ifu, 2i}' for 3 g, then 3  r, then 1 i, 1 ifu, 1 i)

            Notes: the numbers in 'ordering' must be single digit,
                   spec/phot_duration should be duration per exp

        Returns
        -------
        (-1, "ERROR: ...") if there is an issue with the input
        (0, "Request added") if there are no errors
        (0, "Request added, atomicrequests returned ...") if there was an issue with atomicrequest creation
        """
        # TODO: get a better description of cadence/phasesamples/sampletolerance

        # TODO: handle exptime/magnitude in-function?
        requests = self.execute_sql("SELECT object_id, program_id FROM request WHERE status != 'EXPIRED';")
        # check program_id, issue warning if it is a repeat, but allow
        if (pardic['object_id'], pardic['program_id']) in requests:
            # TODO: set a warnings.warn?
            print "program %s has already requested object: %s" % (pardic['program_id'], pardic['object_id'])
            # TODO: require interaction to continue?
        if pardic['object_id'] not in [obj[0] for obj in self.execute_sql('SELECT id FROM object;')]:
            return (-1, "ERROR: object does not exist!")
        # TODO: set default inidate/enddate?

        if 'ordering' in pardic.keys():
            nexpo = pardic['ordering'][1:-1].split(',')
            nexposure = [0, 0, 0, 0, 0]
            for entry in nexpo:
                if entry[1:] == 'ifu':
                    nexposure[0] += int(entry[0])
                elif entry[1:] == 'u':
                    nexposure[1] += int(entry[0])
                elif entry[1:] == 'g':
                    nexposure[2] += int(entry[0])
                elif entry[1:] == 'r':
                    nexposure[3] += int(entry[0])
                elif entry[1:] == 'i':
                    nexposure[4] += int(entry[0])

            # make sure that 'nexposures' and 'ordering' are consistent
            if 'nexposures' in pardic.keys():
                if not '{%s, %s, %s, %s, %s}' % tuple(nexposure) == pardic['nexposures']:
                    return (-1, "ERROR: nexposures and ordering are inconsistent!")
            else:
                pardic['nexposures'] = '{%s, %s, %s, %s, %s}' % tuple(nexposure)

        elif not ('nexposures' in pardic.keys() or 'ordering' in pardic.keys()):
            return (-1, "ERROR: nexposures or ordering is required!")

        keys = list(pardic.keys())
        default_params = ['object_id', 'user_id', 'program_id', 'exptime', 'priority',
                          'inidate', 'enddate']
        for param in default_params:  # check that all required values are provided
            if param not in keys:
                return (-1, "ERROR: %s not in dictionary!" % (param,))
            if not pardic[param]:
                return (-1, "ERROR: no value provided for %s!" % (param,))
        for key in reversed(keys):  # remove any invalid keys
            if key not in ['object_id', 'user_id', 'program_id', 'exptime', 'priority',
                           'inidate', 'enddate', 'marshal_id', 'maxairmass', 'cadence',
                           'phasesamples', 'sampletolerance', 'nexposures', 'ordering']:
                keys.remove(key)
        sql = generate_insert_sql(pardic, keys, 'request')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_request sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_request sql command failed with a ProgrammingError!")
        return (0, "Request added")
        # TODO: test ...

    def update_request(self, pardic):
        """
        Updates the request table with the parameters from the dictionary.

        Parameters
        ----------
        pardic: dict
            required: 'id'
            optional: 'status', 'maxairmass', 'priority', 'inidate', 'enddate'
            Note: 'status' can be 'PENDING', 'ACTIVE', 'COMPLETED', 'CANCELED', or 'EXPIRED'

        Returns
        -------
        (-1, "ERROR: ...") if there was an issue with the updating
        (0, "Requests updated") if the update was successful
        (0, "Requests and atomicrequests updated") if atomicrequests were also updated
        """
        # TODO: if exptime is allowed, significant changes are needed to the atomicrequest update
        # TODO: determine which parameters shouldn't be changed
        keys = list(pardic.keys())
        if 'id' not in keys:
            return (-1, "ERROR: no id provided!")
        if pardic['id'] not in [x[0] for x in self.execute_sql('SELECT id FROM request;')]:
            return (-1, "ERROR: request does not exist!")
        if 'status' in keys:
            if pardic['status'] not in ['PENDING', 'ACTIVE', 'COMPLETED', 'CANCELED', 'EXPIRED']:
                keys.remove('status')
        for key in reversed(keys):  # remove any keys that are invalid or not allowed to be updated
            if key not in ['status', 'maxairmass', 'priority', 'inidate', 'enddate', 'exptime']:
                keys.remove(key)
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        sql = generate_update_sql(pardic, keys, 'request', True)
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_request sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_request sql command failed with a ProgrammingError!")
        # update associated atomicrequests
        update_keys = list(pardic.keys())
        for key in reversed(update_keys):
            if key not in ['priority', 'inidate', 'enddate']:
                keys.remove(key)
        if len(update_keys) == 0:
            return (0, "Requests updated")
        update_sql = "UPDATE atomicrequest SET "
        for param in update_keys:
            if pardic[param]:  # it may be a key with nothing in it
                update_sql += "%s = '%s', " % (param, pardic[param])
        update_sql += "WHERE request_id = %s;" % (pardic['id'],)
        try:
            self.execute_sql(update_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_request atomicrequest sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_request atomicrequest sql command failed with a ProgrammingError!")

        return (0, "Requests and atomicrequests updated")

    def select_from_requests(self, values, where_dict):
        """
        select values from `objects`

        Parameters
        ----------
        values: list of str
            values to be returned
        where_dict: dict
            'param':'value' to be used as WHERE clauses
        values/keys options:
            ['id', 'object_id', 'user_id', 'program_id', 'exptime', 'status', 'priority',
             'inidate', 'enddate', 'marshal_id', 'maxairmass', 'cadence', 'phasesamples',
             'sampletolerance', 'nexposures', 'ordering', 'creationdate', 'lastmodified']

        Returns
        -------
        list: tuples containing the values for each user matching the criteria, empty if no results
        """
        # TODO: test, reconsider return styles
        allowed_params = ['id', 'object_id', 'user_id', 'program_id', 'exptime', 'status', 'priority', 'inidate',
                          'enddate', 'marshal_id', 'maxairmass', 'cadence', 'phasesamples', 'sampletolerance',
                          'nexposures', 'ordering', 'creationdate', 'lastmodified']
        sql = generate_select_sql(values, where_dict, allowed_params, 'request')

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def get_active_requests(self):
        """
        Returns active requests
        Returns
        -------
        list : (object_id, user_id, program_id, marshal_id, exptime, maxairmass, priority, cadence, phasesamples,
             sampletolerance, filters, nexposures, ordering) for each request
        """
        # TODO: test?
        # TODO: move to logic layer
        sql = ("SELECT object_id, user_id, program_id, marshal_id, exptime, maxairmass, priority, cadence, "
               "phasesamples, sampletolerance, filters, nexposures, ordering FROM request WHERE status='ACTIVE';")
        active_requests = self.execute_sql(sql)
        return active_requests

    def expire_requests(self):
        """
        Updates the request table. For all the active requests that were not completed,
            and had an expiry date before than NOW(), are marked as "EXPIRED".
        """
        # TODO: move to logic layer? (requires sql "knowledge")
        # tests written
        sql = "UPDATE request SET status='EXPIRED' WHERE enddate < 'NOW()' AND status != 'COMPLETED';"
        self.execute_sql(sql)
        # TODO: test if the following works
        atomic_sql = ("UPDATE atomicrequest SET status='EXPIRED' WHERE EXISTS (SELECT id FROM request WHERE "
                      "atomicrequest.request_id=request.id AND request.status='EXPIRED')")
        self.execute_sql(atomic_sql)
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

    def add_atomic_request(self, pardic):
        """
        Adds an atomicrequest

        Parameters
        ----------
        pardic: dict
            required:
                'request_id', 'exptime' (duration based on magnitude), 'filter', 'priority', 'inidate', 'enddate'
            optional:
                'object_id', 'order_id' (index of observation order for the request e.g. 0)
            filter options:
                'u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b'

        Returns
        -------
        (-1, "ERROR: ...") if there is an issue
        (0, "Request added") if it succeeded
        """
        # TODO: better description of 'order_id'
        # TODO: determine whether this should handle filter modifications to exptime?
        # TODO: test
        keys = list(pardic.keys())
        for key in ['request_id', 'exptime', 'filter', 'priority', 'inidate', 'enddate']:
            if key not in keys:
                return (01, "ERROR: %s not provided!" % (key,))
            elif not pardic[key]:
                return(-1, "ERROR: no value provided for %s!" % (key,))
        req_obj_stat = self.execute_sql("SELECT object_id, status FROM request "
                                        "WHERE id='%s'" % (pardic['request_id']))
        if not req_obj_stat:  # if there is no request with the id given
            return (-1, "ERROR: request does not exist!")
        if 'object_id' not in keys:
            pardic['object_id'] = req_obj_stat[0][0]
        else:
            if not pardic['object_id'] == req_obj_stat[0][0]:  # check for mismatch of given object_id and request_id
                return (-1, "ERROR: object_id given doesn't match request_id!")
        if req_obj_stat[0][1] == 'EXPIRED':  # check if the request has expired
            return (-1, "ERROR: request has expired!")
        if pardic['filter'] not in ['u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b']:  # check the filter is valid
            return (-1, "ERROR: invalid filter given!")
        for key in reversed(keys):  # remove any invalid keys
            if key not in ['request_id', 'exptime', 'filter', 'priority',
                           'inidate', 'enddate', 'object_id', 'order_id']:
                keys.remove(key)
        sql = generate_insert_sql(pardic, keys, 'atomicrequest')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_atomic_request sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_atomic_request sql command failed with a ProgrammingError!")
        return (0, "Request added")

    def update_atomic_request(self, pardic):
        """
        Updates an atomic request with the parameters from the dictionary

        Parameters
        ----------
        pardic : dict
            required: 'id'
            optional:
                'status', 'priority', 'inidate', 'enddate', 'exptime'
            NOTE: 'status' can be 'PENDING', 'OBSERVED', 'REDUCED', 'EXPIRED' or 'CANCELED'

        Returns
        -------
        (-1, "ERROR: ...") if there was an issue
        (0, "Atomic request updated") if it completed successfully
        """
        # TODO: test, determine which parameters are allowed to be changed
        keys = list(pardic.keys())
        if 'id' not in keys:
            return (-1, "ERROR: no id provided!")
        if 'status' in keys:
            if pardic['status'] not in ['PENDING', 'OBSERVED', 'REDUCED', 'EXPIRED', 'CANCELED']:
                keys.remove('status')  # TODO: remove it, return a -1, or print/warn?
        for key in reversed(keys):  # remove 'id' and any disallowed/invalid keys
            if key not in ['status', 'priority', 'inidate', 'enddate', 'exptime']:
                keys.remove(key)
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")

        sql = generate_update_sql(pardic, keys, 'atomicrequest', lastmodified=True)
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_atomic_request sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_atomic_request sql command failed with a ProgrammingError!")
        return (0, "Atomic request updated")

    def get_request_atomicrequests(self, request_id):
        """
        return the atomicreqests associated with a single request

        Parameters
        ----------
        request_id: int
            The request_id of the desired request

        Returns
        -------
        (id, object_id, order_id, exptime, filter, status, priority )
                for all atomicreqests associated with the desired request
        """
        # TODO: test
        atomic_requests = self.execute_sql("SELECT id, object_id, order_id, exptime, filter, status, priority "
                                           "FROM atomicrequest WHERE request_id = %s" % (request_id,))
        return atomic_requests

    def add_observation_fits(self, header_dict, tel_stats):
        """

        Parameters
        ----------
        header_dict:
        tel_stats:

        Returns
        -------

        """
        if 'atomicrequest_id' in header_dict.keys():
            if not self.execute_sql("SELECT * FROM atomicrequest WHERE id='%s'" % (header_dict['atomicrequest_id'],)):
                return (-1, "ERROR: atomicrequest does not exist!")
            if self.execute_sql("SELECT * FROM observation WHERE atomicrequest_id='%s'"
                                % (header_dict['atomicrequest_id'],)):
                return self.update_observation(header_dict)  # TODO: write update_observation
        atomic_id = header_dict['atomicrequest_id']  # TODO: add ATOM_ID to the header
        # request_id, object_id = self.execute_sql("SELECT request_id, object_id "
        #                                          "FROM atomicrequest WHERE id='%s'" % (atomic_id,))[0]
        # header_dict['request_id'] = request_id
        # header_dict['object_id'] = object_id
        if 'request_id' not in header_dict:  # for testing purposes
            request_id = 1
        if 'object_id' not in header_dict:
            object_id = 1
        # TODO: uncomment above and remove the "request_id=1 and object_id=1" once implemented
        header_keys = list(header_dict.keys())
        for key in ['object_id', 'request_id', 'atomicrequest_id', 'mjd', 'airmass', 'exptime', 'fitsfile', 'lst',
                    'ra', 'dec', 'tel_ra', 'tel_dec', 'tel_az', 'tel_el', 'tel_pa', 'ra_off', 'dec_off']:
            if key not in header_keys:
                return (01, "ERROR: %s not provided!" % (key,))
            elif not header_dict[key]:
                return (-1, "ERROR: no value provided for %s!" % (key,))
        for key in reversed(header_keys):
            if key not in ['object_id', 'request_id', 'atomicrequest_id', 'mjd', 'airmass', 'exptime',
                           'fitsfile', 'imtype', 'lst', 'ra', 'dec', 'tel_ra', 'tel_dec', 'tel_az',
                           'tel_el', 'tel_pa', 'ra_off', 'dec_off', 'camera']:
                header_keys.remove(key)
        sql = generate_insert_sql(header_dict, header_keys, 'observation')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: adding observation sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: adding observation sql command failed with a ProgrammingError!")

        # TODO: see if OUTPUT is able to be used instead of the following
        observation_id = (self.execute_sql("SELECT id FROM observation WHERE atomicrequest_id = '%s'"
                          % (header_dict['atomicrequest_id'],)))[0][0]

        stat_keys = list(tel_stats.keys())
        for key in reversed(stat_keys):
            if key not in ['date', 'dome_status', 'in_temp', 'in_humidity', 'in_dew', 'out_temp', 'out_humidity',
                           'out_dew', 'wind_dir', 'wsp_cur', 'wsp_avg', 'mir_temp', 'top_air', 'pri_temp', 'sec_temp',
                           'flo_temp', 'bot_temp', 'mid_temp', 'top_temp']:
                header_keys.remove(key)
        tel_stats['observation_id'] = int(observation_id)

        stat_sql = generate_insert_sql(tel_stats, stat_keys, 'telescope_stats')
        try:
            self.execute_sql(stat_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: adding tel_stats sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: adding tel_stats sql command failed with a ProgrammingError!")

        # TODO: add other returns for failure cases?
        return (0, "Observation added")
# TODO: write update_observation() and update_telescope_stats()

    def add_reduced_photometry(self, pardic):
        """
        Creates a new object in the phot table

        Parameters
        ----------
        pardic: dict
            required:
                'observation_id', 'astrometry', 'filter', 'reducedfile', 'sexfile', 'biasfile',
                'maskfile', 'flatfile', 'pipeline', 'marshal_phot_id'

        Returns
        -------
        (-1, "ERROR: ...") if there was an issue
        (0, "Photometry added") if the photometry was added successfully
        (0, "Photometry updated for observation_id ...") if the photometry existed and was updated
        """
        # TODO: test
        keys = list(pardic.keys())
        if 'observation_id' not in keys:
            return (-1, "ERROR: observation_id not provided!")
        phot_id = self.execute_sql("SELECT id FROM phot WHERE observation_id='%s'" % (pardic['observation_id']))
        if phot_id:  # if there is already an entry for that observation, update instead
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['astrometry', 'filter', 'reducedfile', 'sexfile', 'biasfile',
                               'maskfile', 'flatfile', 'pipeline', 'marshal_phot_id']:
                    keys.remove(key)
            pardic['id'] = phot_id[0][0]
            update_sql = generate_insert_sql(pardic, keys, 'spec')
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_reduced_photometry update sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_reduced_photometry update sql command failed with a ProgrammingError!")
            return (0, "Photometry updated for observation_id %s" % (pardic['observation_id'],))

        obs = self.execute_sql("SELECT fitsfile FROM observation WHERE id='%s'" % (pardic['observation_id'],))
        if not obs:
            return (-1, "ERROR: no observation with the observation_id")
        else:
            pass
            # TODO: generate the filter here?

        for key in ['observation_id', 'astrometry', 'filter', 'reducedfile', 'sexfile', 'biasfile',
                    'maskfile', 'flatfile', 'pipeline']:  # include 'marshal_phot_id'?
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))
            elif not pardic[key]:
                return (-1, "ERROR: no value provided for %s!" % (key,))

        for key in reversed(keys):  # remove any invalid keys
            if key not in ['observation_id', 'astrometry', 'filter', 'reducedfile', 'sexfile', 'biasfile',
                           'maskfile', 'flatfile', 'pipeline', 'marshal_phot_id']:
                keys.remove(key)
        sql = generate_insert_sql(pardic, keys, 'phot')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_reduced_photometry sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_reduced_photometry sql command failed with a ProgrammingError!")
        # set the atomicrequest's status to 'REDUCED'
        reduced_sql = ("UPDATE atomicrequest SET status='REDUCED' WHERE EXISTS (SELECT id FROM observation "
                       "WHERE observation.atomicrequest_id = atomicrequest.id AND observation.id = '%s';"
                       % (pardic['observation_id'],))  # TODO: test this monstrosity, otherwise can do 2 queries
        try:
            self.execute_sql(reduced_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_reduced_photometry sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_reduced_photometry sql command failed with a ProgrammingError!")
        return (0, "Photometry added")

    def add_reduced_spectrum(self, pardic):
        """
        Creates a new object in the spec table

        Parameters
        ----------
        pardic: dict
            required:
                'observation_id', 'reducedfile', 'sexfile', 'biasfile', 'flatfile', 'imgset',
                'quality', 'cubefile', 'standardfile', 'skysub

        Returns
        -------
        (-1, "ERROR: ...") if there was an issue
        (0, "Spectrum added")  if the spectrum was added successfully
        (0, "Spectrum updated for observation_id ...") if the spectrum existed and was updated
        """
        # TODO: which parameters are required? test
        # TODO: update schedule table indicating reduction?
        keys = list(pardic.keys())
        if 'observation_id' not in keys:
            return (-1, "ERROR: observation_id not provided!")
        spec_id = self.execute_sql("SELECT id FROM spec WHERE observation_id='%s'" % (pardic['observation_id']))
        if spec_id:  # if there is already an entry for that observation, update instead
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['reducedfile', 'sexfile', 'biasfile', 'flatfile', 'imgset', 'quality',
                               'cubefile', 'standardfile', 'marshal_spec_id', 'skysub']:
                    keys.remove(key)
            pardic['id'] = spec_id[0][0]
            update_sql = generate_insert_sql(pardic, keys, 'spec')
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_reduced_spectrum update sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_reduced_spectrum update sql command failed with a ProgrammingError!")
            return (0, "Spectrum updated for observation_id %s" % (pardic['observation_id'],))

        for key in ['observation_id', 'reducedfile', 'sexfile', 'biasfile', 'flatfile', 'imgset', 'quality', 'cubefile',
                    'standardfile', 'skysub']:
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))
            elif not pardic[key]:
                return (-1, "ERROR: no value provided for %s!" % (key,))

        for key in reversed(keys):
            if key not in ['observation_id', 'reducedfile', 'sexfile', 'biasfile', 'flatfile', 'imgset', 'quality', 'cubefile',
                           'standardfile', 'marshal_spec_id', 'skysub']:
                keys.remove(key)
        sql = generate_insert_sql(pardic, keys, 'spec')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_reduced_spectrum sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_reduced_spectrum sql command failed with a ProgrammingError!")
        # set the atomicrequest's status to 'REDUCED'
        reduced_sql = ("UPDATE atomicrequest SET status='REDUCED' WHERE EXISTS (SELECT id FROM observation "
                       "WHERE observation.atomicrequest_id = atomicrequest.id AND observation.id = '%s';"
                       % (pardic['observation_id'],))  # TODO: test this monstrosity, otherwise can do 2 queries
        try:
            self.execute_sql(reduced_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_reduced_spectrum sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_reduced_spectrum sql command failed with a ProgrammingError!")
        return (0, "Spectrum added")

    def add_metrics_phot(self, pardic):
        """
        Creates an entry in metrics_phot or updates an existing metrics entry

        Parameters
        ----------
        pardic: dict
            required: 'phot_id', 'fwhm', 'background', 'zp', 'zperr', 'ellipticity', 'nsources'

        Returns
        -------
        (-1, "ERROR: ...") if there was an issue
        (0, "Photometry metrics updated for phot_id ...") if it updated existing metrics
        (0, "Photometry metrics added") if the metrics were added successfully
        """
        # TODO: which parameters are required? test
        keys = list(pardic.keys())
        if 'phot_id' not in keys:
            return (-1, "ERROR: phot_id not provided!")
        metric_id = self.execute_sql("SELECT id FROM metrics_phot WHERE phot_id='%s'" % (pardic['phot_id']))
        if metric_id:  # if there is already an entry for that observation, update instead
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['fwhm', 'background', 'zp', 'zperr', 'ellipticity', 'nsources']:
                    keys.remove(key)
            pardic['id'] = metric_id[0][0]
            update_sql = generate_insert_sql(pardic, keys, 'metrics_phot')
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_metrics_phot update sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_metrics_phot update sql command failed with a ProgrammingError!")
            return (0, "Photometry metrics updated for phot_id %s" % (pardic['phot_id'],))

        for key in ['fwhm', 'background', 'zp', 'zperr', 'ellipticity', 'nsources']:  # phot_id already tested
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))
            elif not pardic[key]:
                return (-1, "ERROR: no value provided for %s!" % (key,))

        for key in reversed(keys):
            if key not in ['phot_id', 'fwhm', 'background', 'zp', 'zperr', 'ellipticity', 'nsources']:
                keys.remove(key)
        sql = generate_insert_sql(pardic, keys, 'metrics_phot')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_metrics_phot sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_metrics_phot sql command failed with a ProgrammingError!")
        return (0, "Photometry metrics added")

    def add_metrics_spec(self, pardic):
        """
        Creates a new object in the metrics spec stats table with the parameters specified in the dictionary.
        Only one metric shall exist for each observation. If the reduction exists, an update
        is made.
        """
        # TODO: which parameters are required? test
        keys = list(pardic.keys())
        if 'spec_id' not in keys:
            return (-1, "ERROR: spec_id not provided!")
        metric_id = self.execute_sql("SELECT id FROM metrics_spec WHERE spec_id='%s'" % (pardic['spec_id']))
        if metric_id:  # if there is already an entry for that observation, update instead
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['fwhm', 'background', 'line_fwhm']:
                    keys.remove(key)
            pardic['id'] = metric_id[0][0]
            update_sql = generate_insert_sql(pardic, keys, 'metrics_spec')
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_metrics_spec update sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_metrics_spec update sql command failed with a ProgrammingError!")
            return (0, "Spectrum metrics updated for spec_id %s" % (pardic['spec_id'],))

        for key in ['fwhm', 'background', 'line_fwhm']:  # phot_id already tested
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))
            elif not pardic[key]:
                return (-1, "ERROR: no value provided for %s!" % (key,))

        for key in reversed(keys):
            if key not in ['fwhm', 'background', 'line_fwhm']:
                keys.remove(key)
        sql = generate_insert_sql(pardic, keys, 'metrics_spec')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_metrics_spec sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_metrics_spec sql command failed with a ProgrammingError!")
        return (0, "Spectrum metrics added")

    def add_flexure(self, pardic):
        """
        Creates a new object in the flexure table.

        Parameters
        ----------
        pardic: dict
            required: 'rms', 'spec_id_1', 'spec_id_2', 'timestamp1', 'timestamp2'

        Returns
        -------
        (-1, "ERROR: ...") if there was an issue
        (0, "Flexure added")
        """
        # TODO: test
        # TODO: find out what the 'rms' is here
        # TODO: not require timestamp1 and timestamp2, derive them from spec info?
        # if flexure is too high, tell something is wrong?
        keys = list(pardic.keys())
        for key in ['rms', 'spec_id_1', 'spec_id_2', 'timestamp1', 'timestamp2']:
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))
            elif not pardic[key]:
                return (-1, "ERROR: no value provided for %s!" % (key,))
        for key in reversed(keys):
            if key not in ['rms', 'spec_id_1', 'spec_id_2', 'timestamp1', 'timestamp2']:
                keys.remove(key)
        sql = generate_insert_sql(pardic, keys, 'flexure')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_flexure sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_flexure sql command failed with a ProgrammingError!")
        return (0, "Flexure added")

    def add_classification(self, pardic):
        """
        Creates a classification object attached to the reduced spectrum.

        Parameters
        ----------
        pardic: dict
            required: 'spec_id', 'object_id', 'classification', 'redshift', 'redshift_err',  'classifier', 'score'
            optional: 'phase', 'phase_err'

        Returns
        -------
        (-1, "ERROR: ...") if there is an issue
        (0, 'Classification added") if it was successful
        """
        # TODO: clean up the required parameters, test
        keys = list(pardic.keys())
        for key in ['spec_id', 'object_id', 'classification', 'redshift', 'redshift_err', 'classifier', 'score']:
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))
            elif not pardic[key]:
                return (-1, "ERROR: no value provided for %s!" % (key,))
        classified = self.execute_sql("SELECT classification, redshift, redshift_err FROM classification "
                                      "WHERE spec_id = '%s' AND classifier = '%s';" % (pardic['spec_id'],
                                                                                       pardic['classification']))
        if classified:
            return (-1, "ERROR: entry exists for that spectra and classifier with classification %s, redshift %s, "
                        "redshift_err %s. Use `update_classification` if necessary." % (classified[0], classified[1],
                                                                                        classified[2]))
        for key in reversed(keys):  # remove any invalid keys
            if key not in ['spec_id', 'object_id', 'classification', 'redshift', 'redshift_err',  'classifier', 'score',
                           'phase', 'phase_err']:
                keys.remove(key)
        sql = generate_insert_sql(pardic, keys, 'classification')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_classification sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_classification sql command failed with a ProgrammingError!")
        return (0, "Classification added")

    def update_classification(self, pardic):
        """
        Update a classification

        Parameters
        ----------
        pardic: dict
            required: 'id' OR ('spec_id', 'classifier')
            optional:
                'classification', 'redshift', 'redshift_err', 'phase', 'phase_err', 'score'
            Note: this function will not modify id, spec_id, classifier or object_id

        Returns
        -------
        (-1, "ERROR: ...") if there is an issue
        (0, "Classification updated") if it was successful
        """
        # TODO: test, determine which parameters should be allowed
        keys = list(pardic.keys())
        if 'id' in keys:
            try:
                spec_id, classifier = self.execute("SELECT spec_id, classifier FROM classification WHERE id='%s';"
                                                   % (pardic['id']))[0]
            except IndexError:
                return (-1, "ERROR: no classification entry with the given id")
            if 'classifier' in keys:
                if not pardic['classifier'] == classifier:
                    return (-1, "ERROR: classifier provided does not match classification id")
            if 'spec_id' in keys:
                if not pardic['spec_id'] == spec_id:
                    return (-1, "ERROR: spec_id provided does not match classification id")
        elif 'spec_id' in keys and 'classifier' in keys:
            try:
                id = self.execute_sql("SELECT id FROM classification WHERE spec_id='%s' AND "
                                      "classifier='%s'" % (pardic['spec_id'], pardic['classifier']))[0]
            except IndexError:
                return (-1, "ERROR: no classification entry with the given spec_id and classifier")
            pardic['id'] = id[0][0]
        else:
            return (-1, "ERROR: needs id or both spec_id and classifier")

        for key in reversed(keys):  # remove 'id', 'object_id', 'classifier' and any invalid keys
            if key not in ['classification', 'redshift', 'redshift_err', 'phase', 'phase_err', 'score']:
                keys.remove(key)

        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        sql = generate_update_sql(pardic, keys, 'classification')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_classification sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_classification sql command failed with a ProgrammingError!")
        return (0, "Classification updated")


def generate_select_sql(values, where_dict, allowed_params, table):
    """
    generate the sql for a select query

    Parameters
    ----------
    values: list of str
        list of values to return
    where_dict: dict
        {'param':'value',...} adds WHERE param='value'...
    allowed_params: list
        the parameters of the table to be queried
    table: str
        name of the table

    Returns
    -------
    str (sql select query)
    """
    for value in values:
        if value not in allowed_params:
            values.remove(value)
    where_keys = list(where_dict.keys())
    for param in where_keys:
        if param not in allowed_params:
            where_keys.remove(param)
    values_str = ''
    for value in values:
        values_str += value+', '
    values_str = values_str[:-2]

    sql = "SELECT %s FROM %s" % (values_str, table)
    if where_keys:
        sql += "WHERE"
        for key in where_keys:
            sql += " %s = '%s' AND" % (key, where_dict[key])
        sql = sql [:-4]
    return sql


def generate_insert_sql(pardic, param_list, table):
    """
    generate the sql for an insert command

    Parameters
    ----------
    pardic: dict (same as given to upper function)
    param_list: list (list of parameters to insert)
    table: string (name of table)

    Returns
    -------
    sql string
    """
    # TODO: test, re-write insert functions
    columns = "("
    values = "("
    for param in param_list:
        if pardic[param]:  # make sure that there is a value
            columns += (param + ", ")
            values += "'%s', " % (pardic[param],)
    columns = columns[:-2] + ')'
    values = values[:-2] + ')'
    sql = "INSERT INTO %s %s VALUES %s;" % (table, columns, values)
    return sql


def generate_update_sql(pardic, param_list, table, lastmodified=False):
    """
    generate the sql for an update command

    Parameters
    ----------
    pardic: dict (same as given to upper function) (must contain 'id')
    param_list: list (list of parameters to update)
    table: string (name of table)
    lastmodified: bool (if the table has a lastmodified column)

    Returns
    -------
    sql string
    """
    # TODO: test, re-write update functions
    sql = "UPDATE %s SET " % (table,)
    for param in param_list:
        if pardic[param]:  # it may be a key with nothing in it
            sql += "%s = '%s', " % (param, pardic[param])
    if lastmodified:
        sql += "lastmodified = 'NOW()' "
    sql += "WHERE id = %s;" % (pardic['id'],)
    return sql


