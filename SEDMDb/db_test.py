import SedmDb
import SedmDb_tools as db_tools
import numpy as np


db = SedmDb.SedmDB()

# TODO: make tests that cause IntegrityError (giving a value in of the wrong format e.g. string for a decimal parameter) and ProgrammingError (attempting to insert/update a column that doesn't exist)
# TODO: make the functions reject pardic keys that aren't columns

user_dict = {'username': 'test_user', 'name': 'nmo', 'email': 'dodo'}

# there are groups 'default group' and 'second default group'
# (user_id = 1, username=default_user) can be tied to the requests

def test_add_group():  # add_to_group is tested below
    assert db.add_group({'designator': 'default group'}) == (-1, "ERROR: group exists!")
    assert db.add_group({'id': 500}) == (-1, "ERROR: no group designator provided!")


# to test add/remove_user
def test_user_manipulation():
    usernames = [user[0] for user in db.execute_sql('SELECT username from users;')]
    if 'test_user' in usernames:
        db.remove_user(user_dict)    
    added = db.add_user(user_dict)
    # check that the function returned a positive and that it was successful
    assert added[0] == 0
    assert [user[0] for user in db.execute_sql("SELECT username FROM users WHERE username='test_user';")]
    user_id = db.execute_sql("SELECT id FROM users WHERE username='test_user';")[0][0]
    db.add_to_group(user_id, 1)  # add the user to the default group
    # check that it was properly added to its group
    assert (user_id, 1) in db.execute_sql("SELECT user_id, group_id FROM usergroups")
    # check other aspects of adding to groups
    assert db.add_to_group(user_id, 1) == (-1, "ERROR: user already in group!")
    assert db.add_to_group(user_id, 0) == (-1, "ERROR: group does not exist!")
    assert db.add_to_group(0, 1) == (-1, "ERROR: user does not exist!")
    db.add_to_group(user_id, 2)  # multiple groups per user
    assert (user_id, 2) in db.execute_sql("SELECT user_id, group_id FROM usergroups")
    # check that it can't create another user with the same username
    add_again = db.add_user(user_dict)
    assert add_again[0] == -1

    removed = db.remove_user(user_dict)
    # check that the function returned a positive and that it was successful
    assert removed[0] == 0
    assert 'nmo' not in [user[0] for user in db.execute_sql('SELECT name from users WHERE id=2;')]
    # check that the user is removed from all groups
    assert user_id not in [user[0] for user in db.execute_sql("SELECT user_id FROM usergroups")]
    remove_again = db.remove_user(user_dict)
    assert remove_again[0] == -1

    assert db.remove_user({'username': 'unused_name'}) == (-1, "ERROR: no user with that username!")
    assert db.remove_user({'id': 0}) == (-1, "ERROR: no user with that id!")
    assert db.remove_user({'name': 'neo', 'email': '3po'}) == (-1, "ERROR: username or id required!")
    assert db.remove_user({}) == (-1, "ERROR: username or id required!")
    
# TODO: more robust testing of the objects
fixed_object = {'marshal_id': 60, 'name': 'test_obj', 'ra': 20., 'dec': 40., 'typedesig': 'f', 'epoch': 2000.}

def test_object_creation():
    # TODO: test all types of objects
    db.add_object(fixed_object)


request_dict = {'object_id': 2, 'user_id': 1, 'program_id': 1, 'marshal_id': 90, 
                'exptime': '{2700, 600}', 'maxairmass': 3.5, 'priority': 3.5, 'inidate': '2016-11-10',
                'enddate': '2016-11-14', 'nexposures': '{1,2,1,2,0}', 'ordering': '{2g,2r,1ifu,2g,2r}'}
# TODO: add to request_dict to test empty and invalid keys

def test_request_manipulation():
    assert db.add_request(request_dict) == (-1, "ERROR: nexposures and ordering are inconsistent!")
    request_dict.pop('nexposures', None)
    db.add_request(request_dict)  # this re-adds nexposures by generating from ordering
    assert db.execute_sql("SELECT id FROM request WHERE nexposures = '{1,0,4,4,0}' AND status != 'EXPIRED'")
    request_dict.pop('ordering', None)
    request_dict.pop('nexposures', None)
    assert db.add_request(request_dict) == (-1, "ERROR: nexposures or ordering is required!")
    request_dict['ordering'] = '{2g,2r,1ifu,2g,2r}'  # return ordering to the dict

    request_dict.pop('priority', None)
    assert db.add_request(request_dict) == (-1, "ERROR: priority not in dictionary!")
    request_dict['priority'] = 3.5  # re-add it to the dict

    request_dict['object_id'] = 0
    assert db.add_request(request_dict) == (-1, "ERROR: object does not exist!")
    request_dict['object_id'] = 1


def test_request_update():
    obj, status, priority = db.execute_sql("SELECT object_id, status, priority FROM request WHERE id='1'")[0]
    db.update_request({'id': 1, 'object_id': int(obj+1), 'status': 'NEW', 'priority': float(priority+1)})
    new_obj, new_stat, new_priority = db.execute_sql("SELECT object_id, status, priority FROM request WHERE id='1'")[0]
    # object_id isn't allowed to be updated and the given status isn't valid
    assert new_obj == obj and new_stat == status and new_priority == priority+1
    assert db.update_request({'object_id': obj+1, 'status': 'ACTIVE'}) == (-1, "ERROR: no id provided!")
    assert db.update_request({'id': 0, 'maxairmass': 3}) == (-1, "ERROR: request does not exist!")
    assert db.update_request({'id': 1, 'creationdate': '2016-02-04'}) == (-1, "ERROR: no parameters given to update!")
    db.update_request({'id': 1, 'object_id': 1, 'status': 'PENDING', 'priority': 6.})  # reset to defaults


def test_request_expiration():
    test_request = db.execute_sql("SELECT id FROM request WHERE enddate = '2016-11-11';")
    if not test_request:
        request_dict['enddate'] = '2016-11-11'
        db.add_request(request_dict)
    else:
        db.update_request({'id': int(test_request[0][0]), 'status': 'PENDING'})
    assert 'PENDING' == db.execute_sql("SELECT status FROM request WHERE enddate = '2016-11-11';")[0][0]
    db.expire_requests()
    assert 'EXPIRED' == db.execute_sql("SELECT status FROM request WHERE enddate = '2016-11-11';")[0][0]


def test_cancel_request():
    test_request = db.execute_sql("SELECT * FROM request")
    if not test_request:
        db.add_request(request_dict)
    else:
        db.update_request({'id': 1, 'status': 'PENDING'})
    assert 'PENDING' == db.execute_sql("SELECT status FROM request WHERE id=1")[0][0]
    db.cancel_scheduled_request(1)
    assert 'CANCELED' == db.execute_sql("SELECT status FROM request WHERE id=1")[0][0]
    assert db.cancel_scheduled_request(0) == (-1, "ERROR: request does not exist!")

atomic_dict = {'exptime': 1800.0, 'filter': 'f', 'priority': 32., 'inidate': '2016-12-12',
               'enddate': '2017-1-12', 'object_id': 1, 'order_id': 1, 'when': '54'}
# 'when' should do nothing (if it does something sql will error)
# lack of 'request_id' should cause a failure


def test_atomicrequest_manipulation():
    # TODO: test for add
    assert db.add_atomicrequest(atomic_dict) == (-1, "ERROR: request_id not provided!")
    atomic_dict['request_id'] = 3.2
    print db.add_atomicrequest(atomic_dict)
    assert db.add_atomicrequest(atomic_dict) == (-1, "ERROR: request_id must be of type 'int'!")
    atomic_dict['request_id'] = db.get_from_atomicrequests(['request_id'], {'object_id': 1})
    assert db.add_atomicrequest(atomic_dict)[0] == 0


    # TODO: test for update
    pass


def test_request_atomic_requests():
    # create_request_atomic_requests is called by add_request already, so its tests above test valid cases
    assert db_tools.create_request_atomic_requests(0) == (-1, "ERROR: request does not exist!")
    # TODO: find a way to test other situations that doesn't run into (atomicrequests for that request already exist)?


def test_get_request_atomic_requests():
    pass


def test_fits_header_parse():
    fitsfile = '/scr2/sedm/phot/20161012/rc20161012_12_36_34.fits'
    db_tools.add_observation_fitsfile(fitsfile)


def test_add_phot_and_metrics():
    # test add_reduced_photometry   

    # test add_metrics_phot 
    pass


def test_add_spec_and_metrics():
    # test add_reduced_spectrum

    # test add_metrics_spec
    pass


def test_add_flexure():
    pass


def test_classification():
    # test add_classification
    # test adding multiple classifications for one observation (different clssifiers)

    # test update_classification
    pass


if __name__ == '__main__':
    test_add_group()
    test_user_manipulation()
    test_object_creation()
    test_request_manipulation()
    test_request_update()
    test_request_expiration()
    test_atomicrequest_manipulation()

