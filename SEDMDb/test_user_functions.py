import SedmDb

# default entries (shouldn't be changed by functions)
# table: users; 'id': 1, 'username': 'default_user', 'name': 'apo', 'email': 'kiwi'
# table: gropus; 'id': 1, 'designator': 'default group'
# table: groups; 'id': 2, 'designator': 'second default group'
db = SedmDb.SedmDB()


# test add_user
# test remove_user
# test get_from_users


# test add_to_group
# test remove_from_group


# test add_group
def test_add_group():  # add_to_group is tested below
    assert db.add_group({'designator': 'default group'}) == (-1, "ERROR: group exists!")
    assert db.add_group({'id': 500}) == (-1, "ERROR: no group designator provided!")


# to test add/remove_user
def test_user_manipulation():
    user_dict = {'username': 'test_user', 'name': 'nmo', 'email': 'dodo'}
    # test a successful add_user
    added = db.add_user(user_dict)
        # check that the function returned a positive and that it was successful
    assert added[0] == 0
    assert [user[0] for user in db.execute_sql("SELECT username FROM users WHERE username='test_user';")]
    # test get_from_users (no WHERE, WHERE, WHERE with no results)
    assert (db.execute_sql('SELECT username, id, name FROM users') ==
            db.get_from_users(['username', 'id', 'name'], {}))
    assert (db.execute_sql("SELECT username, id, name FROM users WHERE name='nmo'") ==
            db.get_from_users(['username', 'id', 'name'], {'name': 'nmo'}))
    assert (db.execute_sql("SELECT username, id, name FROM users WHERE username='aaa'") ==
            db.get_from_users(['username', 'id', 'name'], {'username': 'aaa'}))
    # test a successful remove_user
    removed = db.remove_user({'username': 'test_user'})
    assert removed[0] == 0
    assert not db.get_from_users(['username'], {'username': 'test_user'})
    # test unsuccessful add_user
    add_again = db.add_user({'username': 'default_user', 'name': 'apo', 'email': 'kiwi'})
    assert add_again[0] == -1


    # test unsuccessful remove_user
    assert db.remove_user({'username': 'unused_name'}) == (-1, "ERROR: no user with that username!")
    assert db.remove_user({'id': 0}) == (-1, "ERROR: no user with that id!")
    assert db.remove_user({'name': 'neo', 'email': '3po'}) == (-1, "ERROR: username or id required!")
    assert db.remove_user({}) == (-1, "ERROR: username or id required!")


    # test for adding a user to a group









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

    removed = db.remove_user(user_dict)
    # check that the function returned a positive and that it was successful
    assert removed[0] == 0
    assert 'nmo' not in [user[0] for user in db.execute_sql('SELECT name from users WHERE id=2;')]
    # check that the user is removed from all groups
    assert user_id not in [user[0] for user in db.execute_sql("SELECT user_id FROM usergroups")]
    remove_again = db.remove_user(user_dict)
    assert remove_again[0] == -1


