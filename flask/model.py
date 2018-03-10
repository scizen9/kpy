import numpy as np
import json
import datetime 
import os, glob
from SEDMDb.SedmDb import SedmDB 
from pandas import DataFrame
from sqlalchemy.exc import IntegrityError
import requests


db = SedmDB(host='localhost', dbname='sedmdb')

def jsonify(dbout, names):
    """
    Receive a response from an SQL DB and converts it into a json file with the names provided.
    """

    if len(dbout) ==0:
        return {}
    elif len(dbout[0]) != len(names):
        return {"message": "Lenght of SQL response different from length of names."}
    else:
        json = []
        for o in dbout:
            j = {}
            for i, n in enumerate(names):
                j[n] = o[i]
            json.append(j)

        return json
            
def updateFollowupConfig(entryUpdate):
    
    message = ""
    if not "Followup" in entryUpdate or not entryUpdate['Followup']:
        message += "STATUS: 422 Invalid Input\n\nFollowup option not specified"
        return
            

    message += "Content-type: text/html"
    message += "\n"
    message += "<title>Test CGI</title>"

    output = open('/home/sedm/kpy/flask/static/test_output%s.txt' % datetime.datetime.utcnow().strftime("%H_%M_%S"),'w')
    
    allFields = ""
    data = json.dumps(entryUpdate)
    
    #for k,v in entryUpdate:
    #    allFields += k + ":" + v + "\n "
    #    output.write(allFields)
    output.write(data)
    output.close()

    return message


def search_stats_file(mydate = None):
    '''
    Returns the last stats file that is present in the system according to the present date.
    It also returns a message stating what date that was.
    '''
    #If the date is specified, we will try to located the right file.
    #None will be returned if it does not exist.
    if ( not mydate is None):
        s= os.path.join("/scr2/sedm/phot", mydate, "stats/stats.log")
        if os.path.isfile(s) and os.path.getsize(s) > 0:
            return s, mydate
        else:
            return None, None

    else:         
        curdate = datetime.datetime.utcnow()
        #Try to find the stat files up to 100 days before today's date.
        i = 0
        while i < 100:
            newdate = curdate - datetime.timedelta(i)
            newdatedir = "%d%02d%02d"%(newdate.year, newdate.month, newdate.day)
            s = os.path.join("/scr2/sedm/phot", newdatedir, "stats/stats.log")
            if os.path.isfile(s) and os.path.getsize(s) > 0:
                return s, newdatedir
            i = i+1
        return None, None


def search_redux_files(mydate = None, user=None):
    '''
    Returns the files that are present in the disk at a given date.
    It also returns a message stating what date that was.

    TODO: This routine that looks in the disk, should look at the database and retrieve only the files which
    correspond to the privileges of the user.
    '''
    #If the date is specified, we will try to located the right file.
    #None will be returned if it does not exist.
    if ( not mydate is None):
        s = os.path.join("/scr2/sedmdrp/redux/", mydate, "*SEDM.txt")
        s2 = os.path.join("/scr2/sedmdrp/redux/", mydate, "image*png")
        s3 = os.path.join("/scr2/sedmdrp/redux/", mydate, "*_SEDM*png")
        s4 = os.path.join("/scr2/sedmdrp/redux/", mydate, "cube*png")
        
        files = glob.glob(s) + glob.glob(s2) + glob.glob(s3) + glob.glob(s4)

        if len(files) == 0:
            files = None

    else:         
        curdate = datetime.datetime.utcnow()
        #Try to find the stat files up to 100 days before today's date.
        i = 0
        files = None
        while i < 100:
            newdate = curdate - datetime.timedelta(i)
            newdatedir = "%d%02d%02d"%(newdate.year, newdate.month, newdate.day)
            s = os.path.join("/scr2/sedmdrp/redux", newdatedir, "*SEDM.txt")
            s2= os.path.join("/scr2/sedmdrp/redux/", newdatedir, "image*png")
            s3 = os.path.join("/scr2/sedmdrp/redux/", newdatedir, "*_SEDM*png")
            s4 = os.path.join("/scr2/sedmdrp/redux/", newdatedir, "cube*png")

            files = glob.glob(s) + glob.glob(s2) + glob.glob(s3) + glob.glob(s4)
            if len(files) > 0:
                mydate = newdatedir
                break

            i = i+1

    if not files is None:
        filenames = [os.path.basename(f) for f in files]
        d = {'filename':filenames}
        df = DataFrame.from_records(d)
    else:
        df = None

    return df, mydate

def search_phot_files(mydate = None, user=None):
    '''
    Returns the files that are present in the disk at a given date.
    It also returns a message stating what date that was.

    TODO: This routine that looks in the disk, should look at the database and retrieve only the files which
    correspond to the privileges of the user.
    '''
    #If the date is specified, we will try to located the right file.
    #None will be returned if it does not exist.
    if ( not mydate is None):
        files = glob.glob(os.path.join("/scr2/sedm/phot/", mydate, "reduced/png/*.png"))
        
        if len(files) == 0:
            files = []

    else:         
        curdate = datetime.datetime.utcnow()
        #Try to find the stat files up to 100 days before today's date.
        i = 0
        files = []
        while i < 100:
            newdate = curdate - datetime.timedelta(i)
            newdatedir = "%d%02d%02d"%(newdate.year, newdate.month, newdate.day)
            files = glob.glob(os.path.join("/scr2/sedm/phot", newdatedir, "reduced/png/*.png"))

            if len(files) > 0:
                mydate = newdatedir
                break

            i = i+1

    files = [os.path.basename(f) for f in files]

    return files, mydate

def get_requests_for_user(user_id):
    '''
    Obtains the DataFrame for the requests that were made:
        - By any member of the group where the user with user_id belongs to.
        - In the last 5 days
    '''
    request_query = ("""SELECT a.designator, o.name, r.inidate, r.enddate, r.priority, r.status, r.lastmodified 
                        FROM request r, object o, allocation a 
                        WHERE o.id = r.object_id AND a.id = r.allocation_id  AND r.enddate > (NOW() - INTERVAL '5 day') AND r.allocation_id IN
                           (SELECT a.id
                            FROM allocation a, groups g, usergroups ug, users u, program p
                            WHERE ug.user_id = u.id AND ug.group_id = g.id AND u.id = %d AND p.group_id = g.id AND a.program_id = p.id
                            ) ORDER BY r.lastmodified DESC;"""% (user_id))
    req = db.execute_sql(request_query)
    req = DataFrame(req, columns=['allocation', 'object', 'start date', 'end date', 'priority','status', 'lastmodified'])

    return req

def get_info_user(username):
    '''
    Returns a json dictionary with the values for the user.
    The fields are: 
        - message: what message we can get from the DB.
        - username:
        - name:
        - email:
        - user_id
        - list of active allocations of the user
        - list of current groups the user belongs to
        - list of groups the user does not belong to
    '''
    if (username==""):
        user_info = None
    elif ('"' in username):
        username = username.split('"')[1]
        user_info = db.execute_sql("SELECT username, name, email, id FROM users WHERE username ='{0}'".format(username))
    else:
        user_info = db.execute_sql("SELECT username, name, email, id FROM users WHERE username LIKE '%{0}%' OR name LIKE '%{1}%' OR email LIKE '%{2}%@%'".format(username, username, username))

    if (not user_info is None and len(user_info)==1):
        user_info = user_info[0]
        username = user_info[0]
        name = user_info[1]
        email = user_info[2]
        user_id = user_info[3]

        user_allocations = db.execute_sql("""SELECT g.designator, p.designator, a.designator, a.inidate, a.enddate  
                                               FROM allocation a, program p, groups g, users u, usergroups ug 
                                               WHERE ug.user_id = u.id AND ug.group_id = g.id 
                                               AND a.program_id=p.id AND p.group_id = g.id AND a.active=TRUE 
                                               AND u.id = %d order by g.designator; """%user_id)
        allocations = jsonify(user_allocations, names=["group", "program", "allocation", "inidate", "enddate"])
        old_groups = db.execute_sql("""SELECT DISTINCT g.designator 
                                        FROM groups g, usergroups ug, users u 
                                        WHERE ug.user_id = u.id AND ug.group_id=g.id AND ug.user_id = u.id AND u.id=%d 
                                        ORDER BY g.designator;"""%user_id)
        new_groups = db.execute_sql("""SELECT DISTINCT g.designator 
                                        FROM groups g  
                                        WHERE g.id NOT IN 
                                            (SELECT DISTINCT group_id 
                                             FROM usergroups ug 
                                             WHERE  ug.user_id=%d) ORDER BY g.designator;"""%user_id)

        u = { "message":"User found.", "username":username, "name":name, "email":email, "id":user_id, \
            "allocations":allocations, "old_groups":old_groups, "new_groups":new_groups}
    elif user_info is None or len(user_info)==0:
        u = {"message": "Found no user with your search criteria. Try with a different name / username / email username."}
    elif len( user_info) > 1:
        users = [u[0] for u in user_info]
        u = {"message": "Found %d users with same search criteria. Use \"\" for exact search. Choose one from: %s"%(len(user_info), users)}
    else:
        u = {"message": "Other issue when searching for the user. %d user found"%len(user_info)}
    return u

def add_group(user_id, g_descriptor):
    '''
    Adds the association between the user and the group in the table usergrouops.
    '''
    group_id= db.execute_sql("SELECT id FROM groups WHERE designator='{0}';".format(g_descriptor))[0][0]
    db.execute_sql("INSERT INTO usergroups values (%d, %d);"%(user_id, group_id))


def remove_group(user_id, g_descriptor):
    '''
        Removes the association between the user and the group in the table usergrouops.
    '''
    group_id= db.execute_sql("SELECT id FROM groups WHERE designator='{0}';".format(g_descriptor))[0][0]
    db.execute_sql("DELETE FROM usergroups ug WHERE ug.user_id=%d AND ug.group_id=%d;"%(user_id, group_id))


def get_p18obsdata(obsdate):
    """
    :param obsdate: Must be in "Year-Month-Day" or "YYYYMMDD" format
    :return: List of dates and average seeing
    """
    #1. Create the URL to get the seeing for the requested night
    p18date = []
    p18seeing = []


    if "-" in obsdate:
        f = datetime.datetime.strptime(obsdate, "%Y-%m-%d") - datetime.timedelta(days=1)
    else:
        f = datetime.datetime.strptime(obsdate, "%Y%m%d") - datetime.timedelta(days=1)

    y, m, d = [f.strftime("%Y"), int(f.strftime("%m")), int(f.strftime("%d"))]
    p18obsdate = "%s-%s-%s" % (y, m, d)

    #2. Get the data from the link
    page = requests.get('http://nera.palomar.caltech.edu/P18_seeing/seeing_log_%s.log' % p18obsdate)
    data = page.content

    #3. Split the page by newlines    
    data = data.split('\n')


    #4. Loop through the data and only use points that have 4 or more seeing values to average
    for i in data:
        i = i.split()

        if len(i) > 5 and int(i[5]) > 4:
            d ='%s %s' %(i[1], i[0])
            p18date.append(datetime.datetime.strptime(d, "%m/%d/%Y %H:%M:%S")
                           + datetime.timedelta(hours=8))
            p18seeing.append(float(i[4]))

    return p18date, p18seeing

def get_allocations_user(user_id):

    res = db.execute_sql(""" SELECT a.designator, p.designator, g.designator, a.time_allocated, a.time_spent
                            FROM allocation a, program p, groups g, usergroups ug
                            WHERE a.program_id = p.id AND p.group_id = g.id 
                            AND g.id = ug.group_id AND a.active is True AND ug.user_id = %d"""%(user_id))

    # create the dataframe and set the allocation names to be linked
    data = DataFrame(res, columns=['allocation', 'program', 'group', 'time allocated', 'time spent'])

    return data


def get_allocations():

    allocations = db.execute_sql("""SELECT a.id, a.designator, p.designator, p.name, a.inidate, a.enddate, a.time_allocated, a.time_spent, a.active
                                               FROM allocation a, program p
                                               WHERE a.program_id=p.id 
                                               ORDER BY a.id DESC; """)
    allocations = jsonify(allocations, names=["id", "name", "program", "description", "inidate", "enddate", "time_allocated", "time_spent", "active"])


    return allocations


def get_programs():

    programs = db.execute_sql("""SELECT p.id, p.designator, p.name, g.designator, p.pi, p.priority
                                               FROM program p, groups g
                                               WHERE p.group_id = g.id 
                                               ORDER BY p.id DESC; """)
    programs = jsonify(programs, names=["id", "designator", "name", "group", "pi", "priority"])


    return programs

def get_all_programs():

    programs = db.get_from_program(["id", "designator"])
    programs = jsonify(programs, names=["id", "name"])


    return programs

def get_all_groups():

    groups = db.execute_sql("SELECT id, designator FROM groups;")
    groups = jsonify(groups, names=["id", "designator"])


    return groups

def delete_allocation(id):

    alloc = db.execute_sql("SELECT * FROM allocation where id=%d;"%id)
    if len(alloc) > 0:
        db.execute_sql("DELETE FROM allocation where id=%d"%id)
        status = 0
        message = "Delected allocation with ID %d"%id
    else:
        status = -1
        message = "No allocation found to delete with ID %d"%id

    return (status, message)


def delete_program(id):

    prog = db.execute_sql("SELECT * FROM program where id=%d;"%id)
    if len(prog) > 0:
        db.execute_sql("DELETE FROM program where id=%d"%id)
        status = 0
        message = "Delected program with ID %d"%id
    else:
        status = -1
        message = "No program found to delete with ID %d"%id

    return (status, message)


def delete_group(id):

    group = db.execute_sql("SELECT * FROM groups where id=%d;"%id)
    if len(group) > 0:
        db.execute_sql("DELETE FROM groups where id=%d"%id)
        status = 0
        message = "Delected group with ID %d"%id
    else:
        status = -1
        message = "No group found to delete with ID %d"%id

    return (status, message)


def get_allocation_stats(user_id):
    """
    Obtains a list of allocations that belong to the user and 
    query the total allocated name and time spent for that allocation.

    If no user_id is provided, all active allocations are returned.
    """
    if (user_id is None):
        res = db.get_from_allocation(["designator", "time_allocated", "time_spent"], {"active":True})
    else:
        res = db.execute_sql(""" SELECT a.designator, a.time_allocated, a.time_spent
                                FROM allocation a, program p, groups g, usergroups ug
                                WHERE a.program_id = p.id AND p.group_id = g.id 
                                AND g.id = ug.group_id AND a.active is True AND ug.user_id = %d"""%(user_id))

    df = DataFrame(res, columns=["designator", "time_allocated", "time_spent"])

    alloc_hours = np.array([ta.total_seconds() / 3600. for ta in df["time_allocated"]])
    spent_hours = np.array([ts.total_seconds() / 3600. for ts in df["time_spent"]])
    free_hours = alloc_hours - spent_hours

    df = df.assign(alloc_hours=alloc_hours, spent_hours=spent_hours, free_hours=free_hours)

    df = df.sort_values(by=["alloc_hours"], ascending=False)

    alloc_names = df["designator"].values
    category = ["alloc_hours", "spent_hours", "free_hours"]


    data = {'allocations' : alloc_names}

    for cat in category:
        data[cat] = df[cat]

    return data

