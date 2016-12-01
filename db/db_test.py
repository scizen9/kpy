import SedmDb
import numpy as np


db = SedmDb.SedmDB()


user_dict = {'group_id': 1, 'username': 'test_user', 'name': 'nmo', 'email': 'dodo'}

if not db.execute_sql("SELECT id FROM groups WHERE id='2'"):  # make sure there are two groups to use
    db.add_group({'designator': 'default group'})
    db.add_group({'designator': 'second default group'})
db.add_user({'group_id': 1, 'username': 'default_user', 'name': 'apo', 'email': 'kiwi'})  # create a user to tie to the requests
default_user_id = db.execute_sql("SELECT id FROM users WHERE username='default_user'")[0][0]


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
fixed_object = {'marshal_id': 90, 'name': 'test_obj', 'ra': 20, 'dec': 40, 'typedesig': 'f', 'epoch': 2000}

if not db.execute_sql("SELECT * FROM object WHERE id='1'"):  # make sure there is a object
    db.add_object(fixed_object)
    fixed_object['marshal_id'] = 60

def test_object_creation():
    # TODO: test all types of objects
    db.add_object(fixed_object, {})


request_dict = {'object_id': 1, 'user_id': default_user_id, 'program_id': 1, 'marshal_id': 90, 
                'exptime': '{2700, 600}', 'maxairmass': '3.5', 'priority': 3.5, 'inidate': '2016-11-10',
                'enddate': '2016-11-14', 'nexposures': '{1,2,1,2,0}', 'ordering': '{2g,2r,1ifu,2g,2r}'}


def test_request_manipulation():
    assert db.add_request(request_dict) == (-1, "ERROR: nexposures and ordering are inconsistent!")
    request_dict.pop('nexposures', None)  # '{1,0,4,3,0}'
    db.add_request(request_dict)
    assert db.execute_sql("SELECT id FROM request WHERE nexposures = '{1,0,4,4,0}' AND status != 'EXPIRED'")
    request_dict.pop('ordering', None)
    assert db.add_request(request_dict) == (-1, "ERROR: nexposures or ordering is required!")
    request_dict['ordering'] = '{2g,2r,1ifu,2g,2r}'  # return ordering to the dict

    request_dict.pop('priority', None)
    assert db.add_request(request_dict) == (-1, "ERROR: priority not in dictionary!")
    request_dict['priority'] = 3.5  # re-add it to the dict

    request_dict['object_id'] = 0
    assert db.add_request(request_dict) == (-1, "ERROR: object does not exist") ????
    request_dict['object_id'] = 1


# TODO: add tests of update_request etc.
def test_request_update():
    object, status, priority = db.execute_sql("SELECT object_id, status, priority FROM request WHERE id='1'")

    pass


def test_request_expiration(request):
    test_request = db.execute_sql("SELECT id FROM request WHERE enddate = '2016-11-11';")
    if not test_request:
        request['enddate'] = '2016-11-11'
        db.add_request(request)
    else:
        db.update_request({'id': test_request[0][0], 'status': 'PENDING'})
    assert 'PENDING' == db.execute_sql("SELECT status FROM request WHERE enddate = '2016-11-11'")[0][0]
    db.expire_requests()
    assert 'EXPIRED' == db.execute_sql("SELECT status FROM request WHERE enddate = '2016-11-11'")[0][0]

def test_fits_header_parse():
    fitsfile= '/scr2/sedm/phot/20161012/rc20161012_12_36_34.fits'
    db.add_observation_fits(fitsfile)

if __name__ == '__main__':
    test_add_group()
    test_user_manipulation()
    test_object_creation()
    test_request_manipulation()
    test_request_expiration(request_dict)

