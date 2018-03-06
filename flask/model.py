import numpy as np
import json
from datetime import datetime, timedelta
import os, glob
from SEDMDb.SedmDb import SedmDB 
from pandas import DataFrame
from sqlalchemy.exc import IntegrityError

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
            return s
        else:
            return None

    else:         
        curdate = datetime.utcnow()
        #Try to find the stat files up to 100 days before today's date.
        i = 0
        while i < 100:
            newdate = curdate - timedelta(i)
            newdatedir = "%d%02d%02d"%(newdate.year, newdate.month, newdate.day)
            s = os.path.join("/scr2/sedm/phot", newdatedir, "stats/stats.log")
            if os.path.isfile(s) and os.path.getsize(s) > 0:
                return s
            i = i+1
        return None


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
        curdate = datetime.utcnow()
        #Try to find the stat files up to 100 days before today's date.
        i = 0
        files = None
        while i < 100:
            newdate = curdate - timedelta(i)
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
        curdate = datetime.utcnow()
        #Try to find the stat files up to 100 days before today's date.
        i = 0
        files = []
        while i < 100:
            newdate = curdate - timedelta(i)
            newdatedir = "%d%02d%02d"%(newdate.year, newdate.month, newdate.day)
            files = glob.glob(os.path.join("/scr2/sedm/phot", newdatedir, "reduced/png/*.png"))

            if len(files) > 0:
                mydate = newdatedir
                break

            i = i+1

    files = [os.path.basename(f) for f in files]

    return files, mydate


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


