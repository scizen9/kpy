import SedmDb
import numpy as np


db = SedmDb.SedmDB()


user_dict = {}
user_dict['group_id'] = 1
user_dict['name'] = 'nmo' 
user_dict['email'] = 'dodo'

#db.add_group({'designator': 'testing group'}) # create default group
#db.add_user({'group_id': 1, 'name': 'apo', 'email': 'kiwi'})  # create a user to tie to the requests

# to test add/remove_user
# db.add_user(user_dict)
# db.remove_user(user_dict)


def test_user_manipulation():
    users = [user[0] for user in db.execute_sql('SELECT name from users;')]
    if 'nmo' in users:
        db.remove_user(user_dict)    
    added = db.add_user(user_dict)    
    assert added[0] == 0
    assert [user[0] for user in db.execute_sql("SELECT name FROM users WHERE name='nmo';")]
    add_again = db.add_user(user_dict)    
    assert add_again[0] == -1

    removed = db.remove_user(user_dict)
    assert removed[0] == 0
    assert 'nmo' not in [user[0] for user in db.execute_sql('SELECT name from users WHERE id=2;')]
    remove_again = db.remove_user(user_dict)
    assert remove_again[0] == -1


    
    # TODO: test for valid/invalid inputs. create an update_user_email?
    
object_dict = {}
object_dict['id'] = 0
object_dict['marshal_id'] = 90
object_dict['name'] = 'test_obj'
object_dict['ra'] = 20
object_dict['dec'] = 40
object_dict['typedesig'] = 'f'
object_dict['epoch'] = 2000

def test_object_creation():
    # TODO: test all types of objects
    db.add_object(object_dict, {})
    

request_dict= {'id': 2, 'object_id': 1, 'user_id': 5, 'program_id': 1, 'marshal_id': 90, 
               'exptime': 2700, 'maxairmass': '3.5', 'status': '', 'priority': 3.5, 'inidate': '2016-11-10',
               'enddate': '2016-12-10', 'cadence': 0, 'phasesamples': '', 'sampletolerance': '',
               'filters': '', 'nexposures': '', 'ordering': '{"2g","1r","1s","2g","2r"}', 'creationdate': '',
               'lastmodified': ''}
def test_request_manipulation():
    db.add_request(request_dict)
    assert db.execute_sql('SELECT * FROM request;')
    # TODO: add tests of update_request etc.



if __name__ == '__main__':
    test_user_manipulation()
    test_object_creation()
    test_request_manipulation()

