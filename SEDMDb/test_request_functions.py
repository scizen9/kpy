import SedmDb

# default entries (shouldn't be changed by functions)
# table: users; 'id': 1, 'username': 'default_user', 'name': 'apo', 'email': 'kiwi'
# table: gropus; 'id': 1, 'designator': 'default group'
# table: groups; 'id': 2, 'designator': 'second default group'
# table: object; 'id': 1, 'name': 'vega', 'ra': 279.23, 'dec': 38.783, 'typedesig': 'f', 'epoch': 2000.
# table: request; 'id': 1, 'object_id': 1, 'user_id': 1, 'exptime': '{ 2400, 360}', 'priority': 6.,
#                 'inidate': '2016-12-30', 'enddate': '2017-02-12', 'nexposures': '{1, 2, 2, 2, 2}'
# table: atomicrequest; 'id': 1 through 'id': 9



# test add_request
# test update_request
# test get_from_requests
# test expire_requests
# test add_atomicrequest
# test update_atomicrequest
# test get_from_atomicrequests