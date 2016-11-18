import SedmDb
import numpy as np


db = SedmDb.SedmDB()


user_dict = {'group_id': 1, 'group_id': 1, 'username': 'test_user', 'name': 'nmo', 'email': 'dodo'}

db.add_group({'designator': 'default group'}) # create default group
db.add_user({'group_id': 1, 'group_id': 1, 'username': 'default_user', 'name': 'apo', 'email': 'kiwi'})  # create a user to tie to the requests
default_user_id = db.execute_sql("SELECT id FROM users WHERE username='default_user'")[0][0]

# to test add/remove_user
def test_user_manipulation():
    users = [user[0] for user in db.execute_sql('SELECT username from users;')]
    if 'test_user' in users:
        db.remove_user(user_dict)    
    added = db.add_user(user_dict)    
    assert added[0] == 0
    assert [user[0] for user in db.execute_sql("SELECT username FROM users WHERE username='test_user';")]
    add_again = db.add_user(user_dict)    
    assert add_again[0] == -1

    removed = db.remove_user(user_dict)
    assert removed[0] == 0
    assert 'nmo' not in [user[0] for user in db.execute_sql('SELECT name from users WHERE id=2;')]
    remove_again = db.remove_user(user_dict)
    assert remove_again[0] == -1

    # TODO: test for valid/invalid inputs. create an update_user?
    
object_dict = {'marshal_id': 90, 'name': 'test_obj', 
               'ra': 20, 'dec': 40, 'typedesig': 'f', 'epoch': 2000}


def test_object_creation():
    # TODO: test all types of objects
    db.add_object(object_dict, {})
    

request_dict = {'object_id': 1, 'user_id': default_user_id, 'program_id': 1, 'marshal_id': 90, 
               'exptime': '{2700, 600}', 'maxairmass': '3.5', 'status': '', 'priority': 3.5, 'inidate': '2016-11-10',
               'enddate': '2016-12-10', 'cadence': 0, 'phasesamples': '', 'sampletolerance': '',
               'nexposures': '{1,1,1,2,0}', 'ordering': '{"2g","1r","1s","2g","2r"}'}


def test_request_manipulation():
    db.add_request(request_dict)    
    assert db.execute_sql('SELECT * FROM request;')
    # TODO: add tests of update_request etc.


def test_request_update():
    pass


def test_request_expiration(request):
    test_request = db.execute_sql("SELECT id FROM request WHERE enddate = '2016-11-11';")
    if not test_request:
        request['enddate'] = '2016-11-11'
        db.add_request(request)
    else:
        db.update_request({'id':test_request[0][0], 'status': 'PENDING'})
    assert 'PENDING' == db.execute_sql("SELECT status FROM request WHERE enddate = '2016-11-11'")[0][0]
    db.expire_requests()
    assert 'EXPIRED' == db.execute_sql("SELECT status FROM request WHERE enddate = '2016-11-11'")[0][0]

def test_fits_header_parse():
    fitsfile= '/scr2/sedm/phot/20161012/rc20161012_12_36_34.fits'
    db.add_observation_fits(fitsfile)

if __name__ == '__main__':
    test_user_manipulation()
    test_object_creation()
    test_request_manipulation()
    test_request_expiration(request_dict)

