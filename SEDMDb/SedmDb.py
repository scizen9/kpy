import sqlalchemy.pool as pool
from sqlalchemy import exc, create_engine
import psycopg2
import numpy as np
import subprocess
import warnings
from astropy.time import Time


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

        Args:
            pardic (dict):
                required:
                    'username' (str),
                    'name' (str),
                    'email'(str)

        Returns:
            (-1, "ERROR...") if there is an issue

            (0, "User added") if the user was added
        """
        # no need to check parameter value types as they are all strings
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
        sql = _generate_insert_sql(pardic, keys, 'users')
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

        Args:
            pardic (dict):
                required:
                    'username' (str)
                    OR
                    'id': (str)

        Returns:
            (-1, "ERROR...") if there was an issue (user doesn't exist, not enough information in pardic)

            (0, "User removed") if the removal was successful
        """
        if 'username' in pardic.keys():
            user_id = self.get_from_users(['id'], {'username': pardic['username']})
            if user_id:
                if user_id[0] == -1:  # if get_from_users failed
                    return user_id
                self.execute_sql("DELETE FROM usergroups WHERE user_id='%s'" % (user_id[0][0],))
                self.execute_sql("DELETE FROM users WHERE id='%s';" % (user_id[0][0],))
                return (0, "User removed")
            else:
                return (-1, "ERROR: no user with that username!")
        elif 'id' in pardic.keys():
            if not (isinstance(pardic['id'], int) or isinstance(pardic['id'], long)):
                return (-1, "ERROR: id must be of type 'int'")
            if pardic['id'] in [x[0] for x in self.execute_sql('SELECT id FROM users;')]:
                self.execute_sql("DELETE FROM usergroups WHERE user_id='%s'" % (pardic['id'],))
                self.execute_sql("DELETE FROM users WHERE id='%s';" % (pardic['id'],))
                return (0, "User removed")
            else:
                return (-1, "ERROR: no user with that id!")
        else:
            return (-1, "ERROR: username or id required!")

    def get_from_users(self, values, where_dict):
        """
        select values from `users`

        Args:
            values: list of str
                values to be returned
            where_dict (dict):
                ``'param':'value'`` to be used as WHERE clauses

            values/keys options:
                'id' (int),
                'username' (str),
                'name' (str),
                'email' (str)]

        Returns:
            list of tuples containing the values for each user matching the criteria
            
            empty list if no users match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """

        # TODO: test, reconsider return styles
        allowed_params = {'id': int, 'username': str, 'name': str, 'email': str}
        sql = _generate_select_sql(values, where_dict, allowed_params, 'users')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)
        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_group(self, pardic):
        """
        Adds a new group. Checks for duplicates in designator

        Args:
            pardic (dict):
                required:
                    'designator' (str)

        Returns:
            (-1, "ERROR...") if no designator was provided or there is already a group with it

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

    def add_usergroup(self, user, group):
        """
        Adds the user as member of the group. Checks for duplicates in name.

        Args:
            user (int):
                id of the user in the 'users' Table
            group: int
                id of the group in the 'groups' Table

        Returns (int):
            (-1, "ERROR...") if there was areason for failure

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

    def remove_from_group(self, user, group):
        """
        removes the user from the group.

        Args:
            user (int):
                id of the user in the 'users' Table
            group: int
                id of the group in the 'groups' Table

        Returns (int):
            (-1, "ERROR...") if there was areason for failure

            (0, "User removed from group") if the removal was successful
        """
        if user not in [user_id[0] for user_id in self.execute_sql('SELECT id FROM users')]:
            return (-1, "ERROR: user does not exist!")
        if group not in [group_id[0] for group_id in self.execute_sql('SELECT id FROM groups')]:
            return (-1, "ERROR: group does not exist!")
        usergroups = self.execute_sql('SELECT user_id, group_id FROM usergroups')
        if (user, group) not in usergroups:
            return (-1, "ERROR: user not in group!")
        else:
            sql = "DELETE FROM usergroups WHERE user_id='%s' AND group_id='%s'" % (user, group)
            try:
                self.execute_sql(sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_to_group sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_to_group sql command failed with a ProgrammingError!")
            return (0, "User removed from group")

    def get_from_usergroups(self, values, where_dict):
        """
        select values from `objects`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options: ['user_id', 'group_id']

        Returns:
            list of tuples containing the values for each user matching the criteria

            empty list if no objects match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
                """
        # TODO: test, reconsider return styles
        allowed_params = {'user_id': int, 'group_id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params, 'usergroups')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_object(self, pardic):
        """
        Creates a new object

        Args:
            pardic (dict):
                required:
                    'name' (str),
                    'typedesig' (str),
                required for a fixed object:
                    'ra' (float) ra in degrees,
                    'dec' (float) dec in degrees,
                    'epoch' (float)
                optional:
                    'iauname' (str),
                    'marshal_id' (int)

                'typedesig' should be one of:
                    'f' (fixed), 'P' (built-in planet or satellite name), 'e' (heliocentric elliptical),
                    'h' (heliocentric hyperbolic), 'p' (heliocentric parabolic), 'E' (geocentric elliptical)

        Returns:
            (-1, "ERROR...") if it failed to add

            (0, "Object added") if the object is added successfully
        """
        param_types = {'name': str, 'typedesig': str, 'ra': float, 'dec': float, 'epoch': float,
                       'iauname': str, 'marshal_id': int}

        obj_keys = list(pardic.keys())
        # TODO: have it update existing object if it already exists?
        if 'marshal_id' in obj_keys:
            if pardic['marshal_id'] in [obj[0] for obj in self.execute_sql('SELECT marshal_id FROM object')]:
                return (-1, "ERROR: object exists!")

        for key in ['name', 'typedesig']:  # check if 'name' and 'typedesig' are provided
            if key not in obj_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(obj_keys):  # remove any extraneous keys
            if key not in ['name', 'typedesig', 'ra', 'dec', 'epoch', 'marshal_id', 'iauname']:
                obj_keys.remove(key)
        type_check = _data_type_check(obj_keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        pardic['name'] = pardic['name'].lower()  # make all of the names the same format for consistant searching
        if pardic['typedesig'] == 'f':
            for key in ['ra', 'dec', 'epoch']:
                if key not in obj_keys:
                    return (-1, "ERROR: %s not provided!" % (key,))
            # dup = self.execute_sql("SELECT id, name FROM object WHERE q3c_radial_query(ra, dec, '%s', '%s', .000278)"
            #                        % (pardic['ra'], pardic['dec']))
            # if dup:  # if there is already an object within an arcsecond
            #     return (-1, "ERROR: there is already an object within 1 arcsec of given coordinates with "
            #                 "id: %s, name: %s" % (object[0][0], object[0][1]))

            obj_sql = _generate_insert_sql(pardic, obj_keys, 'object')
            try:
                self.execute_sql(obj_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_object sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_object sql command failed with a ProgrammingError!")
            return (0, "Fixed object added")
        else:
            obj_sql = _generate_insert_sql(pardic, obj_keys, 'object')
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

    def get_from_object(self, values, where_dict):
        """
        select values from `objects`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options: ['id', 'marshal_id', 'name', 'iauname', 'ra', 'dec', 'typedesig', 'epoch']

        Returns:
            list of tuples containing the values for each user matching the criteria

            empty list if no objects match ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        # TODO: test, reconsider return styles
        allowed_params = {'id': int, 'marshal_id': int, 'name': str, 'iauname': str, 'ra': float, 'dec': float,
                          'typedesig': str, 'epoch': float}
        sql = _generate_select_sql(values, where_dict, allowed_params, 'object')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def get_objects_near(self, ra, dec, radius):
        """
        get object entries with coordinates within a radius

        Args:
            ra (float or int): ra in degrees of the object
            dec (float or int): dec in degrees of the object
            radius (float or int): radius in arcseconds

        Returns:
            list of tuples [(id, name, epoch)] for each object in the radius

            [] if there are no objects found

            (-1, "ERROR...") if there is an issue
        """
        return (-1, "ERROR: q3c is not working")
        
        if not (isinstance(ra, float) or isinstance(ra, int)):
            return (-1, "ERROR: parameter ra must be of type 'float' or type 'int'!")
        if not (isinstance(dec, float) or isinstance(dec, int)):
            return (-1, "ERROR: parameter dec must be of type 'float' or type 'int'!")
        if not (isinstance(radius, float) or isinstance(radius, int)):
            return (-1, "ERROR: parameter radius must be of type 'float' or type 'int'!")

        objects = self.execute_sql("SELECT id, name FROM object WHERE q3c_radial_query(ra, dec, '%s', '%s', '%s')"
                                   % (ra, dec, .000278*radius))
        return objects

    def get_object_id_from_name(self, object_name):
        """
        finds the id of an object given its name or part of its name

        Args:
            object_name (str):

        Returns:
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

    def add_elliptical_heliocentric(self, orbit_params):
        """
        Adds the orbit parameters for an elliptical heliocentric orbit

        Args:
            orbit_params (dict):
                required:
                    'object_id' (int),
                    'inclination' (float),
                    'longascnode_O' (float) (lon. of ascending node),
                    'perihelion_o' (float) (arg. of perihelion),
                    'a' (float) (mean distance AU),
                    'n' (float) (mean daily motion deg/day),
                    'e' (float) (eccentricity),
                    'M' (float) (mean anomaly),
                    'mjdepoch' (int) (epoch, time of 'M'),
                    'D' (int) (equinox year),
                    'M1' (float),
                    'M2' (float) (first and second components of magnitude model)
                optional:
                    's' (float) (angular size at 1 AU)

        Returns:

        """
        param_types = {'object_id': int, 'inclination': float, 'longascnode_O': float, 'perihelion_o': float,
                       'a': float, 'n': float, 'e': float, 'M': float, 'mjdepoch': int, 'D': int, 'M1': float,
                       'M2': float, 's': float}
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['inclination', 'longascnode_O', 'perihelion_o', 'a', 'n', 'e',
                    'M', 'mjdepoch', 'D', 'M1', 'M2', 'object_id']:
            if key not in orb_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(orb_keys):
            if key not in ['inclination', 'longascnode_O', 'perihelion_o', 'a', 'n', 'e',
                           'M', 'mjdepoch', 'D', 'M1', 'M2', 's', 'object_id']:
                orb_keys.remove(key)
        type_check = _data_type_check(orb_keys, orbit_params, param_types)
        if type_check:
            return (-1, type_check)

        orb_sql = _generate_insert_sql(orbit_params, orb_keys, 'elliptical_heliocentric')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_elliptical_orbit sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_elliptical_orbit sql command failed with a ProgrammingError!")
        return (0, "Elliptical heliocentric orbit added")

    def get_from_elliptical_heliocentric(self, values, where_dict):
        """
        select values from `elliptical_heliocentric`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'object_id' (int),
                'inclination' (float),
                'longascnode_O' (float) (lon. of ascending node),
                'perihelion_o' (float) (arg. of perihelion),
                'a' (float) (mean distance AU),
                'n' (float) (mean daily motion deg/day),
                'e' (float) (eccentricity),
                'M' (float) (mean anomaly),
                'mjdepoch' (int) (epoch, time of 'M'),
                'D' (int) (equinox year),
                'M1' (float),
                'M2' (float) (first and second components of magnitude model),
                's' (float) (angular size at 1 AU)

        Returns:
            list of tuples containing the values for each orbit matching the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        # TODO: test, reconsider return styles
        allowed_params = {'object_id': int, 'inclination': float, 'longascnode_O': float, 'perihelion_o': float,
                          'a': float, 'n': float, 'e': float, 'M': float, 'mjdepoch': int, 'D': int, 'M1': float,
                          'M2': float, 's': float, 'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params, 'elliptical_heliocentric')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_hyperbolic_heliocentric(self, orbit_params):
        """
        Adds the orbit parameters for an hyperbolic heliocentric orbit

        Args:
            orbit_params (dict):
                required:
                    'object_id' (int),
                    'T' ('year-month-day' or `~astropy.time.Time` object),
                    'inclination' (float),
                    'longascnode_O' (float) (lon. of ascending node),
                    'perihelion_o' (float) (arg. of perihelion),
                    'e' (float) (eccentricity),
                    'q' (float) (perihelion distance AU),
                    'D' (int) (equinox year),
                    'M1' (float),
                    'M2' (float) (first and second components of magnitude model)

                optional:
                    's' (float) (angular size at 1 AU)

        Returns:

        """
        param_types = {'object_id': int, 'T': 'date', 'e': float, 'inclination': float, 'longascnode_O': float,
                       'perihelion_o': float, 'q': float, 'D': int, 'M1': float, 'M2': float, 's': float}
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['T', 'inclination', 'longascnode_O', 'perihelion_o', 'e', 'q', 'D',
                    'M1', 'M2', 'object_id']:
            if key not in orb_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(orb_keys):
            if key not in ['T', 'inclination', 'longascnode_O', 'perihelion_o', 'e', 'q', 'D',
                           'M1', 'M2', 's', 'object_id']:
                orb_keys.remove(key)
        type_check = _data_type_check(orb_keys, orbit_params, param_types)
        if type_check:
            return (-1, type_check)

        orb_sql = _generate_insert_sql(orbit_params, orb_keys, 'hyperbolic_heliocentric')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_hyperbolic_orbit sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_hyperbolic_orbit sql command failed with a ProgrammingError!")
        return (0, "Hyperbolic heliocentric orbit added")

    def get_from_hyperbolic_heliocentric(self, values, where_dict):
        """
        select values from `hyperbolic_heliocentric`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'object_id' (int),
                'T' ('year-month-day' or `~astropy.time.Time` object),
                'inclination' (float),
                'longascnode_O' (float) (lon. of ascending node),
                'perihelion_o' (float) (arg. of perihelion),
                'e' (float) (eccentricity),
                'q' (float) (perihelion distance AU),
                'D' (int) (equinox year),
                'M1' (float),
                'M2' (float) (first and second components of magnitude model),
                's' (float) (angular size at 1 AU)

        Returns:
            list of tuples containing the values for each orbit matching the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        # TODO: test, reconsider return styles
        allowed_params = {'object_id': int, 'T': 'date', 'e': float, 'inclination': float, 'longascnode_O': float,
                          'perihelion_o': float, 'q': float, 'D': int, 'M1': float, 'M2': float, 's': float, 'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params, 'hyperbolic_heliocentric')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_parabolic_heliocentric(self, orbit_params):
        """
        Adds the orbit parameters for an parabolic heliocentric orbit

        Args:
            orbit_params (dict):
                required:
                    'object_id'(int),
                    'T' (date str),
                    'inclination' (float),
                    'perihelion_o' (float) (arg. of perihelion),
                    'q' (float) (perihelion distance),
                    'longascnode_O' (float) (lon. of ascending node),
                    'D' (int) (equinox year),
                    'M1' (float),
                    'M2' (float) (first and second components of magnitude model)

                optional:
                    's' (float) (angular size at 1 AU)

        Returns:

        """
        param_types = {'object_id': int, 'T': 'date', 'inclination': float, 'longascnode_O': float,
                       'perihelion_o': float, 'q': float, 'D': int, 'M1': float, 'M2': float, 's': float}
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['T', 'inclination', 'longascnode_O', 'perihelion_o', 'q', 'D',
                    'M1', 'M2', 'object_id']:
            if key not in orb_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(orb_keys):
            if key not in ['T', 'inclination', 'longascnode_O', 'perihelion_o', 'q', 'D',
                           'M1', 'M2', 's', 'object_id']:
                orb_keys.remove(key)
        type_check = _data_type_check(orb_keys, orbit_params, param_types)
        if type_check:
            return (-1, type_check)

        orb_sql = _generate_insert_sql(orbit_params, orb_keys, 'parabolic_heliocentric')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_parabolic_orbit sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_parabolic_orbit sql command failed with a ProgrammingError!")
        return (0, "Parabolic heliocentric orbit added")

    def get_from_parabolic_heliocentric(self, values, where_dict):
        """
        select values from `parabolic_heliocentric`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'object_id'(int),
                'T' (date str),
                'inclination' (float),
                'perihelion_o' (float) (arg. of perihelion),
                'q' (float) (perihelion distance),
                'longascnode_O' (float) (lon. of ascending node),
                'D' (int) (equinox year),
                'M1' (float),
                'M2' (float) (first and second components of magnitude model),
                's' (float) (angular size at 1 AU)

        Returns:
            list of tuples containing the values for each orbit matching the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        # TODO: test, reconsider return styles
        allowed_params = {'object_id': int, 'T': 'date', 'inclination': float, 'longascnode_O': float,
                          'perihelion_o': float, 'q': float, 'D': int, 'M1': float, 'M2': float, 's': float, 'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params, 'parabolic_heliocentric')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_earth_satellite(self, orbit_params):
        """
        Adds the orbit parameters for an Earth satellite orbit

        Args:
            orbit_params (dict):
                required:
                    'object_id' (int),
                    'T' ('year-month-day') (epoch of other fields),
                    'inclination' (float),
                    'ra' (float) (ra of ascending node),
                    'e' (float) (eccentricity),
                    'pedigree' (float) (arg. of pedigree),
                    'M' (float) (mean anomaly),
                    'n' (float) (mean motion, revs/day),
                    'decay' (float) (orbit decay rate, rev/day^2),
                    'reforbit' (int) (integral reference orbit number at epoch),

                optional:
                    'drag' (float) (drag coefficient, 1/(Earth radii))

        Returns:

        """
        param_types = {'object_id': int, 'T': 'date', 'e': float, 'inclination': float, 'ra': float,
                       'pedigree': float, 'M': float, 'n': float, 'decay': float, 'reforbit': int, 'drag': float}
        # TODO: query associated table for object already existing, test
        orb_keys = list(orbit_params.keys())
        for key in ['T', 'inclination', 'ra', 'e', 'pedigree', 'M', 'n',
                    'decay', 'reforbit', 'object_id']:
            if key not in orb_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(orb_keys):
            if key not in ['T', 'inclination', 'ra', 'e', 'pedigree', 'M', 'n',
                           'decay', 'reforbit', 'drag', 'object_id']:
                orb_keys.remove(key)
        type_check = _data_type_check(orb_keys, orbit_params, param_types)
        if type_check:
            return (-1, type_check)

        orb_sql = _generate_insert_sql(orbit_params, orb_keys, 'earth_satellite')
        try:
            self.execute_sql(orb_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_earth_satellite_orbit sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_earth_satellite_orbit sql command failed with a ProgrammingError!")
        return (0, "Earth satellite orbit added")

    def get_from_earth_satellite(self, values, where_dict):
        """
        select values from `earth_satellite`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'object_id' (int),
                'T' ('year-month-day') (epoch of other fields),
                'inclination' (float),
                'ra' (float) (ra of ascending node),
                'e' (float) (eccentricity),
                'pedigree' (float) (arg. of pedigree),
                'M' (float) (mean anomaly),
                'n' (float) (mean motion, revs/day),
                'decay' (float) (orbit decay rate, rev/day^2),
                'reforbit' (int) (integral reference orbit number at epoch),
                'drag' (float) (drag coefficient, 1/(Earth radii))

        Returns:
            list of tuples containing the values for each orbit matching the criteria,

            empty list if no orbits match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        # TODO: test, reconsider return styles
        allowed_params = {'object_id': int, 'T': 'date', 'e': float, 'inclination': float, 'ra': float,
                          'pedigree': float, 'M': float, 'n': float, 'decay': float, 'reforbit': int,
                          'drag': float, 'id': int}
        sql = _generate_select_sql(values, where_dict, allowed_params, 'earth_satellite')
        print sql
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def _add_planet_satellite_orbit(self, orbit_params):
        raise NotImplementedError
        # TODO: actually implement? maybe let higher level do this
        # TODO: query ??? for if it already exists

    def add_request(self, pardic):
        """
        Add a request

        Args:
            pardic (dict):
                required:
                    'object_id' (int),
                    'user_id' (int),
                    'program_id' (int),
                    'exptime' (str '{spec_duration, phot_duration}'),
                    'priority' (float),
                    'inidate' ('year-month-day') (start of observing window),
                    'enddate' ('year-month-day') (end of observing window),
                    'nexposures' or 'ordering' (below)

                optional:
                    'marshal_id' (int),
                    'maxairmass' (float) (max allowable airmass for observation, default 2.5),
                    'cadence' (float) (time between periods),
                    'phasesamples' (float) (how many samples in a period),
                    'sampletolerance' (float) (how much tolerance in when the samples can be taken),
                    'nexposures' (str '{# of ifu, # of u, # of g, # of r, # of i}'),
                    'ordering' (str e.g. '{3g, 3r, 1i, 1ifu, 2i}' for 3 of g, then 3 of r, then 1 i, 1 ifu, 1 i)

                Note:
                    the numbers in 'ordering' must be single digit,
                    spec/phot_duration should be duration per exp

        Returns:
            (-1, "ERROR...") if there is an issue with the input

            (0, "Request added") if there are no errors

            (0, "Request added, atomicrequests returned ...") if there was an issue with atomicrequest creation
        """
        # TODO: get a better description of cadence/phasesamples/sampletolerance
        param_types = {'object_id': int, 'user_id': int, 'program_id': int, 'exptime': str, 'priority': float,
                       'inidate': 'date', 'enddate': 'date', 'marshal_id': int, 'maxairmass': float, 'cadence': float,
                       'phasesamples': float, 'sampletolerance': float, 'nexposures': str, 'ordering': str}
        # TODO: handle exptime/magnitude in-function?
        requests = self.execute_sql("SELECT object_id, program_id FROM request WHERE status != 'EXPIRED';")
        # check program_id, issue warning if it is a repeat, but allow
        if (pardic['object_id'], pardic['program_id']) in requests:
            # TODO: set a warnings.warn?
            print "program %s has already requested object: %s" % (pardic['program_id'], pardic['object_id'])
            # TODO: require interaction to continue?
        if pardic['object_id'] not in [obj[0] for obj in self.execute_sql('SELECT id FROM object;')]:
            return (-1, "ERROR: object does not exist!")
        if pardic['user_id'] not in [user[0] for user in self.execute_sql('SELECT id FROM users;')]:
            return (-1, "ERROR: user does not exist!")
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
        for key in reversed(keys):  # remove any invalid keys
            if key not in ['object_id', 'user_id', 'program_id', 'exptime', 'priority',
                           'inidate', 'enddate', 'marshal_id', 'maxairmass', 'cadence',
                           'phasesamples', 'sampletolerance', 'nexposures', 'ordering']:
                keys.remove(key)
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'request')
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

        Args:
            pardic (dict):
                required:
                    'id' (int)
                optional:
                    'status' (str),
                    'maxairmass' (float),
                    'priority' (float),
                    'inidate' ('year-month-day'),
                    'enddate' ('year-month-day')
                Note: 'status' can be 'PENDING', 'ACTIVE', 'COMPLETED', 'CANCELED', or 'EXPIRED'

        Returns:
            (-1, "ERROR...") if there was an issue with the updating

            (0, "Requests updated") if the update was successful

            (0, "Requests and atomicrequests updated") if atomicrequests were also updated
        """
        # TODO: if exptime is allowed, significant changes are needed to the atomicrequest update
        # TODO: determine which parameters shouldn't be changed
        param_types = {'id': int, 'status': str, 'maxairmass': float, 'priority': float,
                       'inidate': 'date', 'enddate': 'date'}
        keys = list(pardic.keys())
        if 'id' not in keys:
            return (-1, "ERROR: no id provided!")
        elif not (isinstance(pardic['id'], int) or isinstance(pardic['id'], long)):
            return (-1, "ERROR: parameter id must be of type 'int'!")
        if pardic['id'] not in [x[0] for x in self.execute_sql('SELECT id FROM request;')]:
            return (-1, "ERROR: request does not exist!")
        if 'status' in keys:
            if pardic['status'] not in ['PENDING', 'ACTIVE', 'COMPLETED', 'CANCELED', 'EXPIRED']:
                keys.remove('status')
        for key in reversed(keys):  # remove any keys that are invalid or not allowed to be updated
            if key not in ['status', 'maxairmass', 'priority', 'inidate', 'enddate']:
                keys.remove(key)
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_update_sql(pardic, keys, 'request', True)
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
                update_keys.remove(key)
        if len(update_keys) == 0:
            return (0, "Requests updated")

        update_sql = _generate_update_sql(pardic, update_keys, 'atomicrequest')
        try:
            self.execute_sql(update_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_request atomicrequest sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_request atomicrequest sql command failed with a ProgrammingError!")

        return (0, "Requests and atomicrequests updated")

    def get_from_request(self, values, where_dict):
        """
        select values from `request`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'object_id' (int),
                'user_id' (int),
                'program_id' (int),
                'exptime' (str),
                'priority' (float),
                'inidate' ('year-month-day'),
                'enddate' ('year-month-day'),
                'marshal_id' (int),
                'maxairmass' (float),
                'cadence' (float),
                'phasesamples' (float),
                'sampletolerance' (float),
                'filters' (str),
                'nexposures' (str),
                'ordering' (str)
                'status' (str),
                'creationdate' ('year-month-day'),
                'lastmodified' ('year-month-day')

        Returns:
            list of tuples containing the values for each request matching the criteria,

            empty list if no requests match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        # TODO: test, reconsider return styles
        allowed_params = {'id': int, 'object_id': int, 'user_id': int, 'program_id': int, 'exptime': str, 'status': str,
                          'priority': float, 'inidate': 'date', 'enddate': 'date', 'marshal_id': int,
                          'maxairmass': float, 'cadence': float, 'phasesamples': float, 'sampletolerance': float,
                          'filters': str, 'nexposures': str, 'ordering': str, 'creationdate': 'date',
                          'lastmodified': 'date'}
        sql = _generate_select_sql(values, where_dict, allowed_params, 'request')
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def expire_requests(self):
        """
        Updates the request table. For all the active requests that were not completed,
            and had an expiry date before than NOW(), are marked as "EXPIRED".

        Returns:
            (0, "Requests expired")
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

    def add_atomicrequest(self, pardic):
        """
        Adds an atomicrequest

        Args:
            pardic (dict):
                required:
                    'request_id' (int),
                    'exptime' (float) (duration based on magnitude/filter),
                    'filter' (str),
                    'priority' (float),
                    'inidate' ('year-month-day'),
                    'enddate' ('year-month-day')
                optional:
                    'object_id' (int),
                    'order_id' (int) (index of observation order for the request e.g. 1)
                filter options:
                    'u', 'g', 'r', 'i', 'ifu', 'ifu_a', 'ifu_b'

        Returns:
            (-1, "ERROR...") if there is an issue

            (0, "Request added") if it succeeded
        """
        # TODO: better description of 'order_id'
        # TODO: determine whether this should handle filter modifications to exptime?
        # TODO: test
        param_types = {'request_id': int, 'exptime': float, 'filter': str, 'priority': float, 'inidate': 'date',
                       'enddate': 'date', 'object_id': int, 'order_id': int}
        keys = list(pardic.keys())
        for key in ['request_id', 'exptime', 'filter', 'priority', 'inidate', 'enddate']:
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))

        req_obj_stat = self.get_from_request(['object_id', 'status'], {'id': pardic['request_id']})
        if not req_obj_stat:  # if there is no request with the id given
            return (-1, "ERROR: request does not exist!")
        elif req_obj_stat[0] == -1:
            return req_obj_stat

        if 'object_id' not in keys:
            pardic['object_id'] = int(req_obj_stat[0][0])
            keys.append('object_id')
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
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)
        sql = _generate_insert_sql(pardic, keys, 'atomicrequest')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_atomic_request sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_atomic_request sql command failed with a ProgrammingError!")
        return (0, "Request added")

    def update_atomicrequest(self, pardic):
        """
        Updates an atomic request with the parameters from the dictionary

        Args:
            pardic (dict):
                required:
                    'id' (int)
                optional:
                    'status' (str),
                    'priority' (float),
                    'inidate' ('year-month-day'),
                    'enddate' ('year-month-day'),
                    'exptime' (float)
                NOTE: 'status' can be 'PENDING', 'OBSERVED', 'REDUCED', 'EXPIRED' or 'CANCELED'

        Returns:
            (-1, "ERROR...") if there was an issue

            (0, "Atomic request updated") if it completed successfully
        """
        param_types = {'id': int, 'exptime': float, 'status': str, 'priority': float, 'inidate': 'date',
                       'enddate': 'date'}
        # TODO: test, determine which parameters are allowed to be changed
        keys = list(pardic.keys())
        if 'id' not in keys:
            return (-1, "ERROR: no id provided!")
        elif pardic['id'] not in [x[0] for x in self.execute_sql('SELECT id FROM atomicrequest;')]:
            return (-1, "ERROR: atomicrequest does not exist!")

        if 'status' in keys:
            if pardic['status'] not in ['PENDING', 'OBSERVED', 'REDUCED', 'EXPIRED', 'CANCELED']:
                keys.remove('status')  # TODO: remove it, return a -1, or print/warn?
        for key in reversed(keys):  # remove 'id' and any disallowed/invalid keys
            if key not in ['status', 'priority', 'inidate', 'enddate', 'exptime']:
                keys.remove(key)
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_update_sql(pardic, keys, 'atomicrequest', lastmodified=True)
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_atomic_request sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_atomic_request sql command failed with a ProgrammingError!")
        return (0, "Atomicrequest updated")

    def get_from_atomicrequest(self, values, where_dict):
        """
        select values from `atomicrequest`

        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'object_id' (int),
                'request_id' (int),
                'order_id' (int),
                'exptime' (str),
                'filter' (str),
                'status' (str),
                'priority' (float),
                'inidate' ('year-month-day'),
                'enddate' ('year-month-day'),
                'creationdate' ('year-month-day'),
                'lastmodified' ('year-month-day')

        Returns:
            list of tuples containing the values for each atomicrequest matching the criteria

            empty list if no atomicrequests match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        # TODO: test, reconsider return styles
        allowed_params = {'id': int, 'object_id': int, 'request_id': int, 'order_id': int, 'exptime': float,
                          'filter': str, 'status': str, 'priority': float, 'inidate': 'date', 'enddate': 'date',
                          'creationdate': 'date', 'lastmodified': 'date'}

        sql = _generate_select_sql(values, where_dict, allowed_params, 'atomicrequest')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_observation(self, header_dict):
        """
        Adds an observation

        Args:
            header_dict (dict):
                required:
                    'object_id' (int),
                    'request_id' (int),
                    'atomicrequest_id' (int),
                    'mjd' (float),
                    'airmass' (float),
                    'exptime' (float),
                    'fitsfile' (str),
                    'lst' (str),
                    'ra' (float),
                    'dec' (float),
                    'tel_ra' (str),
                    'tel_dec' (str),
                    'tel_az' (float),
                    'tel_el' (float),
                    'tel_pa' (float),
                    'ra_off' (float),
                    'dec_off' (float),
                optional:
                    'imtype' (str),
                    'camera' (str)

        Returns:
            (-1, "ERROR...") if there was an issue

            (0, "Observation added") if it completed successfully
        """
        header_types = {'object_id': int, 'request_id': int, 'atomicrequest_id': int, 'mjd': float, 'airmass': float,
                        'exptime': float, 'fitsfile': str, 'lst': str, 'ra': float, 'dec': float, 'tel_ra': str,
                        'tel_dec': str, 'tel_az': float, 'tel_el': float, 'tel_pa': float, 'ra_off': float,
                        'dec_off': float, 'imtype': str, 'camera': str}
        if 'atomicrequest_id' in header_dict.keys():
            if not (isinstance(header_dict['atomicrequest_id'], int) or
                        isinstance(header_dict['atomicrequest_id'], long)):  # prevent the select sql from failing
                return (-1, "ERROR: atomicrequest must be of type 'int'!")
            if not self.execute_sql("SELECT * FROM atomicrequest WHERE id='%s'" % (header_dict['atomicrequest_id'],)):
                return (-1, "ERROR: atomicrequest does not exist!")
            elif self.execute_sql("SELECT * FROM observation WHERE atomicrequest_id='%s'"
                                  % (header_dict['atomicrequest_id'],)):
                return self.update_observation(header_dict)  # TODO: write update_observation
        else:
            return (-1, "ERROR: no atomicrequest_id provided!")

        header_keys = list(header_dict.keys())
        for key in ['object_id', 'request_id', 'atomicrequest_id', 'mjd', 'airmass', 'exptime', 'fitsfile', 'lst',
                    'ra', 'dec', 'tel_ra', 'tel_dec', 'tel_az', 'tel_el', 'tel_pa', 'ra_off', 'dec_off']:
            if key not in header_keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(header_keys):
            if key not in ['object_id', 'request_id', 'atomicrequest_id', 'mjd', 'airmass', 'exptime',
                           'fitsfile', 'imtype', 'lst', 'ra', 'dec', 'tel_ra', 'tel_dec', 'tel_az',
                           'tel_el', 'tel_pa', 'ra_off', 'dec_off', 'camera']:
                header_keys.remove(key)
        type_check = _data_type_check(header_keys, header_dict, header_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(header_dict, header_keys, 'observation')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: adding observation sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: adding observation sql command failed with a ProgrammingError!")

        # TODO: add other returns for failure cases?
        return (0, "Observation added")

    def update_observation(self, pardic):
        """

        Args:
            pardic (dict):
                required:
                    'id' (int)
                    OR
                    'atomicrequest_id' (int)
                optional:
                    'mjd' (float),
                    'airmass' (float),
                    'exptime' (float),
                    'fitsfile' (str),
                    'lst' (str),
                    'ra' (float),
                    'dec' (float),
                    'tel_ra' (str),
                    'tel_dec' (str),
                    'tel_az' (float),
                    'tel_el' (float),
                    'tel_pa' (float),
                    'ra_off' (float),
                    'dec_off' (float),
                    'imtype' (str),
                    'camera' (str)

        Returns:

        """
        # TODO: reconsider allowed parameters
        param_types = {'id': int, 'atomicrequest_id': int, 'mjd': float, 'airmass': float,
                       'exptime': float, 'fitsfile': str, 'lst': str, 'ra': float, 'dec': float, 'tel_ra': str,
                       'tel_dec': str, 'tel_az': float, 'tel_el': float, 'tel_pa': float, 'ra_off': float,
                       'dec_off': float, 'imtype': str, 'camera': str}
        keys = list(pardic.keys())
        if 'id' not in keys and 'atomicrequest_id' not in keys:
            return (-1, "ERROR: neither id nor atomicrequest_id provided!")
        elif 'id' not in keys:
            obs = self.get_from_observation(['id'], {'atomic_request_id': pardic['atomicrequest_id']})
            if not obs:
                return (-1, "ERROR: there is no observation with that atomicrequest_id")
            elif obs[0] == -1:
                return obs
            pardic['id'] = obs[0][0]

        elif pardic['id'] not in [x[0] for x in self.execute_sql('SELECT id FROM observation;')]:
            return (-1, "ERROR: observation does not exist!")

        for key in reversed(keys):  # remove any keys that are invalid or not allowed to be updated
            if key not in ['mjd', 'airmass', 'exptime', 'fitsfile', 'lst', 'ra', 'dec', 'tel_ra',
                           'tel_dec', 'tel_az', 'tel_el', 'tel_pa', 'ra_off', 'dec_off', 'imtype', 'camera']:
                keys.remove(key)
        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_update_sql(pardic, keys, 'observation')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_observation sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_observation sql command failed with a ProgrammingError!")

    def get_from_observation(self, values, where_dict):
        """
        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'object_id' (int),
                'request_id' (int),
                'atomicrequest_id' (int),
                'mjd' (float),
                'airmass' (float),
                'exptime' (float),
                'fitsfile' (str),
                'lst' (str),
                'ra' (float),
                'dec' (float),
                'tel_ra' (str),
                'tel_dec' (str),
                'tel_az' (float),
                'tel_el' (float),
                'tel_pa' (float),
                'ra_off' (float),
                'dec_off' (float),
                'imtype' (str),
                'camera' (str)

        Returns:
            list of tuples containing the values for each observation matching the criteria

            empty list if no observations match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'object_id': int, 'request_id': int, 'atomicrequest_id': int, 'mjd': float, 'airmass': float,
                          'exptime': float, 'fitsfile': str, 'lst': str, 'ra': float, 'dec': float, 'tel_ra': str,
                          'tel_dec': str, 'tel_az': float, 'tel_el': float, 'tel_pa': float, 'ra_off': float,
                          'dec_off': float, 'imtype': str, 'camera': str, 'id': int}

        sql = _generate_select_sql(values, where_dict, allowed_params, 'observation')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_telescope_stats(self, tel_stats):
        """
        Adds telescope stats associated with an observation

        Args:
            tel_stats:
                required:
                    'observation_id' (int),
                    'date' ('year-month-day'),
                    'dome_status' (str),
                    'in_temp' (float),
                    'in_humidity' (float),
                    'in_dew' (float),
                    'out_temp' (float),
                    'out_humidity' (float),
                    'out_dew' (float),
                    'wind_dir' (float),
                    'wsp_cur' (float),
                    'wsp_avg' (float),
                    'mir_temp' (float),
                    'top_air' (float),
                    'pri_temp' (float),
                    'sec_temp' (float),
                    'flo_temp' (float),
                    'bot_temp' (float),
                    'mid_temp' (float),
                    'top_temp' (float)

        Returns:
            (-1, "ERROR...") if an issue occurs

            (0, "Telescope stats added") if successful
        """
        telstat_types = {'date': 'date', 'dome_status': str, 'in_temp': float, 'in_humidity': float, 'in_dew': float,
                         'out_temp': float, 'out_humidity': float, 'out_dew': float, 'wind_dir': float,
                         'wsp_cur': float, 'wsp_avg': float, 'mir_temp': float, 'top_air': float, 'pri_temp': float,
                         'sec_temp': float, 'flo_temp': float, 'bot_temp': float, 'mid_temp': float, 'top_temp': float,
                         'observation_id': int}
        stat_keys = list(tel_stats.keys())
        for key in reversed(stat_keys):
            if key not in ['date', 'dome_status', 'in_temp', 'in_humidity', 'in_dew', 'out_temp', 'out_humidity',
                           'out_dew', 'wind_dir', 'wsp_cur', 'wsp_avg', 'mir_temp', 'top_air', 'pri_temp', 'sec_temp',
                           'flo_temp', 'bot_temp', 'mid_temp', 'top_temp', 'observation_id']:
                stat_keys.remove(key)
        type_check = _data_type_check(stat_keys, tel_stats, telstat_types)
        if type_check:
            return (-1, type_check)

        stat_sql = _generate_insert_sql(tel_stats, stat_keys, 'telescope_stats')
        try:
            self.execute_sql(stat_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: adding tel_stats sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: adding tel_stats sql command failed with a ProgrammingError!")

        return (0, "Telescope stats added")

    # TODO: write update_observation() and update_telescope_stats()

    def get_from_telescope_stats(self, values, where_dict):
        """
        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'observation_id' (int),
                'date' ('year-month-day'),
                'dome_status' (str),
                'in_temp' (float),
                'in_humidity' (float),
                'in_dew' (float),
                'out_temp' (float),
                'out_humidity' (float),
                'out_dew' (float),
                'wind_dir' (float),
                'wsp_cur' (float),
                'wsp_avg' (float),
                'mir_temp' (float),
                'top_air' (float),
                'pri_temp' (float),
                'sec_temp' (float),
                'flo_temp' (float),
                'bot_temp' (float),
                'mid_temp' (float),
                'top_temp' (float)

        Returns:
            list of tuples containing the values for telescope stats matching the criteria

            empty list if no telescope_stats entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'date': 'date', 'dome_status': str, 'in_temp': float, 'in_humidity': float, 'in_dew': float,
                          'out_temp': float, 'out_humidity': float, 'out_dew': float, 'wind_dir': float,
                          'wsp_cur': float, 'wsp_avg': float, 'mir_temp': float, 'top_air': float, 'pri_temp': float,
                          'sec_temp': float, 'flo_temp': float, 'bot_temp': float, 'mid_temp': float, 'top_temp': float,
                          'observation_id': int, 'id': int}

        sql = _generate_select_sql(values, where_dict, allowed_params, 'telescope_stats')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    # TODO: separate add_phot/spec and update_phot/spec
    def add_phot(self, pardic):
        """
        Adds reduced photometry or, if the observation already has photometry, updates it

        Args:
            pardic (dict):
                required:
                    'observation_id' (int),
                    'astrometry' ('true' or 'false'),
                    'filter' (str),
                    'reducedfile' (str),
                    'sexfile' (str),
                    'biasfile' (str),
                    'maskfile' (str),
                    'flatfile' (str),
                    'pipeline' (str),
                    'marshal_phot_id' (int)

        Returns:
            (-1, "ERROR...") if there was an issue

            (0, "Photometry added") if the photometry was added successfully

            (0, "Photometry updated for observation_id ...") if the photometry existed and was updated
        """
        param_types = {'observation_id': int, 'astrometry': 'bool', 'filter': str, 'reducedfile': str, 'sexfile': str,
                       'biasfile': str, 'maskfile': str, 'flatfile': str, 'pipeline': str, 'marshal_phot_id': int}
        # TODO: test
        keys = list(pardic.keys())
        if 'observation_id' not in keys:
            return (-1, "ERROR: observation_id not provided!")
        phot_id = self.get_from_phot(['id'], {'observation_id': pardic['observation_id']})
        if phot_id:  # if there is already an entry for that observation, update instead
            if phot_id[0] == -1:
                return phot_id
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['astrometry', 'filter', 'reducedfile', 'sexfile', 'biasfile',
                               'maskfile', 'flatfile', 'pipeline', 'marshal_phot_id']:
                    keys.remove(key)
            pardic['id'] = phot_id[0][0]
            type_check = _data_type_check(keys, pardic, param_types)
            if type_check:
                return (-1, type_check)

            update_sql = _generate_insert_sql(pardic, keys, 'spec')
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_reduced_photometry update sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_reduced_photometry update sql command failed with a ProgrammingError!")
            return (0, "Photometry updated for observation_id %s" % (pardic['observation_id'],))

        obs = self.get_from_observation(['fitsfile'], {'id': pardic['observation_id']})
        if not obs:
            return (-1, "ERROR: no observation with the observation_id")
        elif obs[0] == -1:
            return obs
        else:
            pass
            # TODO: generate the filter here?

        for key in ['observation_id', 'astrometry', 'filter', 'reducedfile', 'sexfile', 'biasfile',
                    'maskfile', 'flatfile', 'pipeline']:  # include 'marshal_phot_id'?
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))

        for key in reversed(keys):  # remove any invalid keys
            if key not in ['observation_id', 'astrometry', 'filter', 'reducedfile', 'sexfile', 'biasfile',
                           'maskfile', 'flatfile', 'pipeline', 'marshal_phot_id']:
                keys.remove(key)
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'phot')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_reduced_photometry sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_reduced_photometry sql command failed with a ProgrammingError!")
        # set the atomicrequest's status to 'REDUCED'
        reduced_sql = ("UPDATE atomicrequest SET status='REDUCED' WHERE EXISTS (SELECT id FROM observation "
                       "WHERE observation.atomicrequest_id = atomicrequest.id AND observation.id = '%s'"
                       % (pardic['observation_id'],))  # TODO: test this monstrosity, otherwise can do 2 queries
        try:
            self.execute_sql(reduced_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_reduced_photometry sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_reduced_photometry sql command failed with a ProgrammingError!")
        return (0, "Photometry added")

    def get_from_phot(self, values, where_dict):
        """
        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'observation_id' (int),
                'astrometry' ('true' or 'false'),
                'filter' (str),
                'reducedfile' (str),
                'sexfile' (str),
                'biasfile' (str),
                'maskfile' (str),
                'flatfile' (str),
                'pipeline' (str),
                'marshal_phot_id' (int)

        Returns:
            list of tuples containing the values for phot entries matching the criteria

            empty list if no phot entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'observation_id': int, 'astrometry': 'bool', 'filter': str, 'reducedfile': str, 'sexfile': str,
                          'biasfile': str, 'maskfile': str, 'flatfile': str, 'pipeline': str, 'marshal_phot_id': int,
                          'id': int}

        sql = _generate_select_sql(values, where_dict, allowed_params, 'phot')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_spec(self, pardic):
        """
        Adds the reduced spectrum or, if the observation already has a spectrum, updates it

        Args:
            pardic (dict):
                required:
                    'observation_id' (int),
                    'reducedfile' (str),
                    'sexfile' (str),
                    'biasfile' (str),
                    'flatfile' (str),
                    'imgset' (str),
                    'quality' (int),
                    'cubefile' (str),
                    'standardfile' (str),
                    'skysub' ('true' or 'false')

        Returns:
            (-1, "ERROR...") if there was an issue

            (0, "Spectrum added")  if the spectrum was added successfully

            (0, "Spectrum updated for observation_id ...") if the spectrum existed and was updated
        """
        param_types = {'observation_id': int, 'reducedfile': str, 'sexfile': str, 'biasfile': str, 'flatfile': str,
                       'imgset': str, 'quality': int, 'cubefile': str, 'standardfile': str, 'skysub': 'bool'}
        # TODO: which parameters are required? test
        # TODO: update schedule table indicating reduction?
        keys = list(pardic.keys())
        if 'observation_id' not in keys:
            return (-1, "ERROR: observation_id not provided!")
        spec_id = self.get_from_spec(['id'], {'observation_id': pardic['observation_id']})
        if spec_id:  # if there is already an entry for that observation, update instead
            if spec_id[0] == -1:
                return spec_id
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['reducedfile', 'sexfile', 'biasfile', 'flatfile', 'imgset', 'quality',
                               'cubefile', 'standardfile', 'marshal_spec_id', 'skysub']:
                    keys.remove(key)
            pardic['id'] = spec_id[0][0]
            update_sql = _generate_update_sql(pardic, keys, 'spec')
            print update_sql
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_reduced_spectrum update sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_reduced_spectrum update sql command failed with a ProgrammingError!")
            return (0, "Spectrum updated for observation_id %s" % (pardic['observation_id'],))
        obs_id = self.get_from_observation(['id'], {'id': pardic['observation_id']})
        if not obs_id:
            return (-1, "ERROR: no observation exists with the given id!")
        elif obs_id[0] == -1:
            return obs_id

        for key in ['observation_id', 'reducedfile', 'sexfile', 'biasfile', 'flatfile', 'imgset', 'quality', 'cubefile',
                    'standardfile', 'skysub']:
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))

        for key in reversed(keys):
            if key not in ['observation_id', 'reducedfile', 'sexfile', 'biasfile', 'flatfile', 'imgset', 'quality',
                           'cubefile',
                           'standardfile', 'marshal_spec_id', 'skysub']:
                keys.remove(key)
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'spec')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_reduced_spectrum sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_reduced_spectrum sql command failed with a ProgrammingError!")
        # set the atomicrequest's status to 'REDUCED'
        reduced_sql = ("UPDATE atomicrequest SET status='REDUCED' WHERE EXISTS (SELECT id FROM observation "
                       "WHERE observation.atomicrequest_id = atomicrequest.id AND observation.id = '%s');"
                       % (pardic['observation_id'],))  # TODO: test this monstrosity, otherwise can do 2 queries
        try:
            self.execute_sql(reduced_sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_reduced_spectrum sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_reduced_spectrum sql command failed with a ProgrammingError!")
        return (0, "Spectrum added")

    def get_from_spec(self, values, where_dict):
        """
        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'observation_id' (int),
                'reducedfile' (str),
                'sexfile' (str),
                'biasfile' (str),
                'flatfile' (str),
                'imgset' (str),
                'quality' (int),
                'cubefile' (str),
                'standardfile' (str),ear
                'skysub' ('true' or 'false')

        Returns:
            list of tuples containing the values for spectra matching the criteria

            empty list if no spec entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'observation_id': int, 'reducedfile': str, 'sexfile': str, 'biasfile': str, 'flatfile': str,
                          'imgset': str, 'quality': int, 'cubefile': str, 'standardfile': str,
                          'skysub': 'bool', 'id': int}

        sql = _generate_select_sql(values, where_dict, allowed_params, 'spec')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_metrics_phot(self, pardic):
        """
        Creates an entry in metrics_phot or updates an existing metrics entry

        Args:
            pardic (dict):
                required:
                    'phot_id' (int),
                    'fwhm' (float),
                    'background' (float),
                    'zp' (float),
                    'zperr' (float),
                    'ellipticity' (float),
                    'nsources' (int)

        Returns:
            (-1, "ERROR...") if there was an issue

            (0, "Photometry metrics updated for phot_id ...") if it updated existing metrics

            (0, "Photometry metrics added") if the metrics were added successfully
        """
        param_types = {'phot_id': int, 'fwhm': float, 'background': float, 'zp': float,
                       'zperr': float, 'ellipticity': float, 'nsources': int}
        # TODO: which parameters are required? test
        keys = list(pardic.keys())
        if 'phot_id' not in keys:
            return (-1, "ERROR: phot_id not provided!")
        metric_id = self.get_from_metrics_phot(['id'], {'phot_id': pardic['phot_id']})
        if metric_id:  # if there is already an entry for that observation, update instead
            if metric_id[0] == -1:
                return metric_id
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['fwhm', 'background', 'zp', 'zperr', 'ellipticity', 'nsources']:
                    keys.remove(key)
            pardic['id'] = metric_id[0][0]
            type_check = _data_type_check(keys, pardic, param_types)
            if type_check:
                return (-1, type_check)

            update_sql = _generate_insert_sql(pardic, keys, 'metrics_phot')
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_metrics_phot update sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_metrics_phot update sql command failed with a ProgrammingError!")
            return (0, "Photometry metrics updated for phot_id %s" % (pardic['phot_id'],))
        ph_id = self.get_from_phot(['id'], {'id': pardic['phot_id']})
        if not ph_id:
            return (-1, "ERROR: no photometry exists with the given id!")
        elif ph_id[0] == -1:
            return ph_id

        for key in ['fwhm', 'background', 'zp', 'zperr', 'ellipticity', 'nsources']:  # phot_id already tested
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))

        for key in reversed(keys):
            if key not in ['phot_id', 'fwhm', 'background', 'zp', 'zperr', 'ellipticity', 'nsources']:
                keys.remove(key)
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'metrics_phot')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_metrics_phot sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_metrics_phot sql command failed with a ProgrammingError!")
        return (0, "Photometry metrics added")

    def get_from_metrics_phot(self, values, where_dict):
        """
        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'phot_id' (int),
                'fwhm' (float),
                'background' (float),
                'zp' (float),
                'zperr' (float),
                'ellipticity' (float),
                'nsources' (int)

        Returns:
            list of tuples containing the values for metrics matching the criteria

            empty list if no metrics_phot entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'phot_id': int, 'fwhm': float, 'background': float, 'zp': float,
                          'zperr': float, 'ellipticity': float, 'nsources': int, 'id': int}

        sql = _generate_select_sql(values, where_dict, allowed_params, 'metrics_phot')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_metrics_spec(self, pardic):
        """
        Creates a new object in the metrics spec stats table with the parameters specified in the dictionary.
        Only one metric exists for each observation. If the reduction exists, an update is made.

        Args:
            pardic:
                required:
                    'spec_id' (int)
                optional:
                    'fwhm' (float),
                    'background' (float),
                    'line_fwhm' (int)

        Returns:
            (-1: "ERROR...") if there is an issue

            (0, "Spectrum metrics added") if it completes successfully
        """
        param_types = {'spec_id': int, 'fwhm': float, 'background': float, 'line_fwhm': int}
        # TODO: which parameters are required? test
        keys = list(pardic.keys())
        if 'spec_id' not in keys:
            return (-1, "ERROR: spec_id not provided!")
        metric_id = self.get_from_metrics_spec(['id'], {'spec_id': pardic['spec_id']})
        if metric_id:  # if there is already an entry for that observation, update instead
            if metric_id[0] == -1:
                return metric_id
            for key in reversed(keys):  # TODO: test the updating
                if key not in ['fwhm', 'background', 'line_fwhm']:
                    keys.remove(key)
            pardic['id'] = metric_id[0][0]
            type_check = _data_type_check(keys, pardic, param_types)
            if type_check:
                return (-1, type_check)

            update_sql = _generate_insert_sql(pardic, keys, 'metrics_spec')
            try:
                self.execute_sql(update_sql)
            except exc.IntegrityError:
                return (-1, "ERROR: add_metrics_spec update sql command failed with an IntegrityError!")
            except exc.ProgrammingError:
                return (-1, "ERROR: add_metrics_spec update sql command failed with a ProgrammingError!")
            return (0, "Spectrum metrics updated for spec_id %s" % (pardic['spec_id'],))
        sp_id = self.get_from_spec(['id'], {'id': pardic['spec_id']})
        if not sp_id:
            return (-1, "ERROR: no spectrum exists with the given id!")
        elif sp_id[0] == -1:
            return sp_id

        for key in []:  # spec_id already tested
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))

        for key in reversed(keys):
            if key not in ['spec_id', 'fwhm', 'background', 'line_fwhm']:
                keys.remove(key)
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'metrics_spec')
        print sql, type_check
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_metrics_spec sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_metrics_spec sql command failed with a ProgrammingError!")
        return (0, "Spectrum metrics added")

    def get_from_metrics_spec(self, values, where_dict):
        """
        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'spec_id' (int),
                'fwhm' (float),
                'background' (float),
                'line_fwhm' (int)

        Returns:
            list of tuples containing the values for metrics matching the criteria

            empty list if no metrics_spec entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'spec_id': int, 'fwhm': float, 'background': float, 'line_fwhm': int, 'id': int}

        sql = _generate_select_sql(values, where_dict, allowed_params, 'metrics_spec')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_flexure(self, pardic):
        """
        Creates a new object in the flexure table.

        Args:
            pardic (dict):
                required:
                    'rms' (float),
                    'spec_id_1' (int),
                    'spec_id_2' (int),
                    'timestamp1' ('year-month-day hour:minute:second'),
                    'timestamp2' ('year-month-day hour:minute:second')

        Returns:
            (-1, "ERROR...") if there was an issue

            (0, "Flexure added")
        """
        param_types = {'rms': float, 'spec_id_1': int, 'spec_id_2': int,
                       'timestamp1': 'datetime', 'timestamp2': 'datetime'}
        # TODO: test
        # TODO: find out what the 'rms' is here
        # TODO: not require timestamp1 and timestamp2, derive them from spec info?
        # if flexure is too high, tell something is wrong?
        keys = list(pardic.keys())
        for key in ['rms', 'spec_id_1', 'spec_id_2', 'timestamp1', 'timestamp2']:
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        for key in reversed(keys):
            if key not in ['rms', 'spec_id_1', 'spec_id_2', 'timestamp1', 'timestamp2']:
                keys.remove(key)
        sp1_id = self.get_from_spec(['id'], {'id': pardic['spec_id_1']})
        if not sp1_id:
            return (-1, "ERROR: no spectrum exists with the given spec_id_1!")
        elif sp1_id[0] == -1:
            return sp1_id
        sp2_id = self.get_from_spec(['id'], {'id': pardic['spec_id_2']})
        if not sp2_id:
            return (-1, "ERROR: no spectrum exists with the given spec_id_2!")
        elif sp2_id[0] == -1:
            return sp2_id

        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'flexure')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: add_flexure sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: add_flexure sql command failed with a ProgrammingError!")
        return (0, "Flexure added")

    def get_from_flexure(self, values, where_dict):
        """
        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'rms' (float),
                'spec_id_1' (int),
                'spec_id_2' (int),
                'timestamp1' ('year-month-day hour:minute:second'),
                'timestamp2' ('year-month-day hour:minute:second')

        Returns:
            list of tuples containing the values for flexure matching the criteria

            empty list if no flexure entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'id': int, 'rms': float, 'spec_id_1': int, 'spec_id_2': int,
                          'timestamp1': 'datetime', 'timestamp2': 'datetime'}

        sql = _generate_select_sql(values, where_dict, allowed_params, 'flexure')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results

    def add_classification(self, pardic):
        """
        Creates a classification object attached to the reduced spectrum.

        Args:
            pardic (dict):
                required:
                    'spec_id' (int),
                    'object_id' (int),
                    'classification' (str),
                    'redshift' (float),
                    'redshift_err' (float),
                    'classifier' (str),
                    'score' (float)
                optional:
                    'phase' (float),
                    'phase_err' (float)
        Returns:
            (-1, "ERROR...") if there is an issue

            (0, 'Classification added") if it was successful
        """
        param_types = {'spec_id': int, 'object_id': int, 'classification': str, 'redshift': float,
                       'redshift_err': float, 'classifier': str, 'score': float, 'phase': float, 'phase_err': float}
        # TODO: clean up the required parameters, test
        keys = list(pardic.keys())
        for key in ['spec_id', 'object_id', 'classification', 'redshift', 'redshift_err', 'classifier', 'score']:
            if key not in keys:
                return (-1, "ERROR: %s not provided!" % (key,))
        classified = self.get_from_classification(['classification', 'redshift', 'redshift_err'],
                                     {'spec_id': pardic['spec_id'], 'classifier': pardic['classifier']})
        if classified:
            if classified[0] == -1:
                return classified
            return (-1, "ERROR: entry exists for that spectrum and classifier with classification %s, redshift %s, "
                        "redshift_err %s. Use `update_classification` if necessary."
                        % (classified[0][0], classified[0][1], classified[0][2]))
        for key in reversed(keys):  # remove any invalid keys
            if key not in ['spec_id', 'object_id', 'classification', 'redshift', 'redshift_err', 'classifier', 'score',
                           'phase', 'phase_err']:
                keys.remove(key)
        

        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sql = _generate_insert_sql(pardic, keys, 'classification')
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

        Args:
            pardic (dict):
                required:
                    'id' (int)
                     OR
                     'spec_id' (int),
                     'classifier' (str)
                optional:
                    'classification' (str),
                    'redshift' (float),
                    'redshift_err' (float),
                    'phase' (float),
                    'phase_err' (float),
                    'score' (float)
                Note: this function will not modify id, spec_id, classifier or object_id

        Returns:
            (-1, "ERROR...") if there is an issue

            (0, "Classification updated") if it was successful
        """
        param_types = {'id': int, 'spec_id': int, 'classifier': str, 'classification': str, 'redshift': float,
                       'redshift_err': float, 'phase': float, 'phase_err': float, 'score': float}
        # TODO: test, determine which parameters should be allowed
        keys = list(pardic.keys())
        if 'id' in keys:
            id_classifier = self.get_from_classification(['spec_id', 'classifier'], {'id': pardic['id']})
            if not id_classifier:
                return (-1, "ERROR: no classification entry with the given id!")
            elif id_classifier[0] == -1:
                return id_classifier

            if 'classifier' in keys:
                if not pardic['classifier'] == id_classifier[0][1]:
                    return (-1, "ERROR: classifier provided does not match classification id!")
            if 'spec_id' in keys:
                if not pardic['spec_id'] == id_classifier[0][0]:
                    return (-1, "ERROR: spec_id provided does not match classification id!")
        elif 'spec_id' in keys and 'classifier' in keys:
            id = self.get_from_classification(['id'], {'spec_id': pardic['spec_id'],
                                                       'classifier': pardic['classifier']})
            if not id:
                return (-1, "ERROR: no classification entry with the given spec_id and classifier!")
            elif id[0] == -1:
                return id
            pardic['id'] = id[0][0]
        else:
            return (-1, "ERROR: needs id or both spec_id and classifier")

        for key in reversed(keys):  # remove 'id', 'object_id', 'classifier' and any invalid keys
            if key not in ['classification', 'redshift', 'redshift_err', 'phase', 'phase_err', 'score']:
                keys.remove(key)

        if len(keys) == 0:
            return (-1, "ERROR: no parameters given to update!")
        type_check = _data_type_check(keys, pardic, param_types)
        if type_check:
            return (-1, type_check)

        sp_id = self.get_from_spec(['id'], {'id': pardic['spec_id']})
        if not sp_id:
            return (-1, "ERROR: no spectrum exists with the given spec_id!")
        elif sp_id[0] == -1:
            return sp_id
        obj_id = self.get_from_object(['id'], {'id': pardic['object_id']})
        if not obj_id:
            return (-1, "ERROR: no object exists with the given object_id!")
        elif obj_id[0] == -1:
            return obj_id

        sql = _generate_update_sql(pardic, keys, 'classification')
        try:
            self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: update_classification sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: update_classification sql command failed with a ProgrammingError!")
        return (0, "Classification updated")

    def get_from_classification(self, values, where_dict):
        """
        Args:
            values (list): list of str
                values to be returned
            where_dict (dict):
                'param':'value' to be used as WHERE clauses
            values/keys options:
                'id' (int),
                'spec_id' (int),
                'classifier' (str),
                'classification' (str),
                'redshift' (float),
                'redshift_err' (float),
                'phase' (float),
                'phase_err' (float),
                'score' (float)

        Returns:
            list of tuples containing the values for classifications matching the criteria

            empty list if no classification entries match the ``where_dict`` criteria

            (-1, "ERROR...") if there was an issue
        """
        allowed_params = {'spec_id': int, 'object_id': int, 'classification': str, 'redshift': float, 'id': int,
                          'redshift_err': float, 'classifier': str, 'score': float, 'phase': float, 'phase_err': float}

        sql = _generate_select_sql(values, where_dict, allowed_params, 'classification')  # checks type and
        if sql[0] == 'E':  # if the sql generation returned an error
            return (-1, sql)

        try:
            results = self.execute_sql(sql)
        except exc.IntegrityError:
            return (-1, "ERROR: sql command failed with an IntegrityError!")
        except exc.ProgrammingError:
            return (-1, "ERROR: sql command failed with a ProgrammingError!")
        return results


def _data_type_check(keys, pardic, value_types):
    """
    make sure the values given match the data types required in the database

    Args:
        keys (list): keys to be tested
        pardic (dict): keys and keys' values
        value_types (dict): keys and types the values should be (e.g. {'ra': float})

        types are: str, float, int, 'date', 'datetime', 'bool'

    Returns:
        "ERROR..." if a value was of the wrong type

        None if all matched
    """
    for key in keys:
        if value_types[key] == str:
            if pardic[key] is not None:
                pass  # all values are added to the queries as strings anyway
            else:
                return "ERROR: %s must not be None!" % (key,)
        elif value_types[key] == float:
            try:
                pardic[key] = float(pardic[key])
            except ValueError:
                return "ERROR: %s must be of %s!" % (key, str(float)[1:-1])
        elif value_types[key] == int:
            try:
                pardic[key] = int(pardic[key])
            except ValueError:
                return "ERROR: %s must be of %s!" % (key, str(int)[1:-1])
        elif value_types[key] == 'date':
            # TODO: find a better way to check this? does it actually modify pardic for the function?
            try:
                pardic[key] = str(Time(pardic[key])).split(' ')[0]
            except ValueError:
                return "ERROR: %s must be of the format 'year-month-day'!" % (key,)
        elif value_types[key] == 'datetime':
            try:
                pardic[key] = str(Time(pardic[key]))
            except ValueError:
                return "ERROR: %s must be of the format 'year-month-day hour:minute:second'!" % (key,)
        elif value_types[key] == 'bool':
            if not (pardic[key] == 'false' or pardic[key] == 'true'):
                return "ERROR: %s must be either 'true' or 'false'" % (key,)
        elif not isinstance(pardic[key], value_types[key]):
            return "ERROR: parameter %s must be of %s!" % (key, str(value_types[key])[1:-1])
    return None


def _generate_select_sql(values, where_dict, allowed_params, table):
    """
    generate the sql for a select query

    Args:
        values (list): list of names of values
            list of values to return
        where_dict (dict):
            {'param':'value',...} adds WHERE param='value'...
        allowed_params (dict):
            the parameters of the table to be queried and their type {'param':type,...}
        table (str):
            name of the table

    Returns:
        str (sql select query)

        "ERROR..." if they type_check fails
    """
    for value in reversed(values):
        if value not in allowed_params.keys():
            values.remove(value)
    where_keys = list(where_dict.keys())
    for param in reversed(where_keys):
        if param not in allowed_params:
            where_keys.remove(param)
    type_check = _data_type_check(where_keys, where_dict, allowed_params)
    if type_check:
        return type_check
    if not values:
        return "ERROR: no valid values requested!"

    values_str = ''
    for value in values:
        values_str += value + ', '
    values_str = values_str[:-2]

    sql = "SELECT %s FROM %s" % (values_str, table)
    if where_keys:
        sql += " WHERE"
        for key in where_keys:
            sql += " %s = '%s' AND" % (key, where_dict[key])
        sql = sql[:-4] + ";"
    return sql


def _generate_insert_sql(pardic, param_list, table):
    """
    generate the sql for an insert command

    Args:
        pardic (dict): (same as given to calling function)
        param_list (list): list of names of parameters to insert
        table (str): name of table

    Returns:
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


def _generate_update_sql(pardic, param_list, table, lastmodified=False):
    """
    generate the sql for an update command

    Args:
        pardic (dict): (same as given to upper function) (must contain 'id')
        param_list (list): list of names of parameters to update
        table (str): name of table
        lastmodified (bool): if the table has a lastmodified column

    Returns:
        sql string
    """
    # TODO: test, re-write update functions
    sql = "UPDATE %s SET" % (table,)
    for param in param_list:
        if pardic[param]:  # it may be a key with nothing in it
            sql += " %s = '%s'," % (param, pardic[param])
    if lastmodified:
        sql += " lastmodified = 'NOW()'"
    else:
        sql = sql[:-1]
    sql += " WHERE id = %s;" % (pardic['id'],)
    return sql
