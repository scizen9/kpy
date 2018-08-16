import numpy as np
import json
import datetime 
import os, glob
from SEDMDb.SedmDb import SedmDB
from pandas import DataFrame
import psycopg2
from psycopg2 import extras

from sqlalchemy.exc import IntegrityError
import requests
import SEDMrph.fitsutils as fitsutils

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

def get_request_by_id(request_id):
    """
    Grab pertinent information on the request id and create a page that can
    be used to update the request
    :param request_id: 
    :return: 
    """

    request_query = """SELECT r.id as req_id, r.object_id as obj_id, 
                            r.user_id, r.marshal_id, r.exptime, r.maxairmass,
                            r.max_fwhm, r.min_moon_dist, r.max_moon_illum, 
                            r.max_cloud_cover, r.status, 
                            r.priority as reqpriority, r.inidate, r.enddate,
                            r.cadence, r.phasesamples, r.sampletolerance, 
                            r.filters, r.nexposures, r.obs_seq, r.seq_repeats,
                            r.seq_completed, r.last_obs_jd, r.creationdate,
                            r.lastmodified, r.allocation_id, r.marshal_id, 
                            o.id as objid, o.name as objname, o.iauname, o.ra, o."dec",
                            o.typedesig, o.epoch, o.magnitude, o.creationdate, 
                            u.id as user_id, u.email, a.id as all_id, a.inidate, a.enddate, a.time_spent, 
                            a.time_allocated, a.program_id, a.active, 
                            p.designator, p.name, p.group_id, p.pi,
                            p.time_allocated, r.priority, p.inidate,
                            p.enddate 
                            FROM "public".request r
                            INNER JOIN "public"."object" o ON ( r.object_id = o.id  )  
                            INNER JOIN "public".users u ON ( r.user_id = u.id  )  
                            INNER JOIN "public".allocation a ON ( r.allocation_id = a.id  )  
                            INNER JOIN "public".program p ON ( a.program_id = p.id  )
                            WHERE r.id = %s
                            """ % request_id

    conn = psycopg2.connect("dbname=sedmdb user=sedmuser "
                            "host=localhost " 
                            "password=user$edm1235")
    cursor = conn.cursor(cursor_factory=psycopg2.extras.DictCursor)
    cursor.execute(request_query)
    results = cursor.fetchone()

    if isinstance(results, psycopg2.extras.DictRow):
        obs_dict = parse_obs(results['obs_seq'], results['exptime'])
    else:
        obs_dict = ''

    return results, obs_dict

def parse_obs(sequence, exptime_sequence):
    """
    Parse all available filters
    :param seq: 
    :param exptime: 
    :return: 
    """
    flt_list = ['ifu', 'r', 'g', 'i', 'u']
    seq = list(sequence)
    exptime = list(exptime_sequence)

    obs_dict = {}
    print seq, exptime
    for flt in flt_list:
        index = [i for i, s in enumerate(seq) if flt in s]

        if index:
            obs_dict['%s_checked' % flt] = 'checked'
            obs_dict['%s_exptime' % flt] = exptime[index[0]]
            obs_dict['%s_repeat' % flt] = int(seq[index[0]].replace(flt,""))

            for i in index:
                seq.pop(index[0])
                exptime.pop(index[0])

        else:
            obs_dict['%s_checked' % flt] = ''
            obs_dict['%s_exptime' % flt] = 0
            obs_dict['%s_repeat' % flt] = 0
    return obs_dict

def delete_request_by_id(id):
    """
    
    :param id: 
    :return: 
    """
    obs_dict = {'id': int(id),
                'status': 'CANCELED'}

    ret = db.update_request(obs_dict)
    print ret
    return "Canceled"


def parse_update(update_dict):
    """
    
    :param update_dict: 
    :return: 
    """

    if update_dict['type'] == 'priority':
        print int(update_dict['id'])
        print db.update_request({'id': int(update_dict['id']),
                           'priority': update_dict['priority']})
    elif update_dict['type'] == 'status':
        print db.update_request({'id': int(update_dict['id']),
                           'status': update_dict['status']})
    elif update_dict['type'] == 'filters':
        pass

    return {'response': True}

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

        basedir = os.path.join("/scr2/sedmdrp/redux/", mydate)
        patterns = ["*SEDM.txt", "image*png", "*_SEDM*png", "cube*png", "Standard_Correction.png", \
            "spec_*.png", "spec_*.txt", "ifu_spaxels*.png", "*flat3d.png", "*wavesolution_dispersionmap.png"]

        files = []

        for p in patterns:
            files.extend(glob.glob(os.path.join(basedir, p)))

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
            basedir = os.path.join("/scr2/sedmdrp/redux/", newdatedir)

            patterns = ["*SEDM.txt", "image*png", "*_SEDM*png", "cube*png", "Standard_Correction.png", \
                "spec_*.png", "spec_*.txt", "ifu_spaxels*.png", "*flat3d.png", "*wavesolution_dispersionmap.png"]

            files = []

            for p in patterns:
                files.extend(glob.glob(os.path.join(basedir, p)))

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
        filesraw = glob.glob(os.path.join("/scr2/sedm/phot/", mydate, "pngraw/*all.png"))

        files = files + filesraw

        if len(files) == 0:
            files = None

    else:         
        curdate = datetime.datetime.utcnow()
        #Try to find the stat files up to 100 days before today's date.
        i = 0
        files = []
        while i < 100:
            newdate = curdate - datetime.timedelta(i)
            newdatedir = "%d%02d%02d"%(newdate.year, newdate.month, newdate.day)
            files = glob.glob(os.path.join("/scr2/sedm/phot", newdatedir, "reduced/png/*.png"))
            filesraw = glob.glob(os.path.join("/scr2/sedm/phot", newdatedir, "pngraw/*all.png"))
            files = files + filesraw

            if len(files) > 0:
                mydate = newdatedir
                break

            i = i+1

    if not files is None:
        files.sort(reverse=True)
        d = {'filename':files}
        df = DataFrame.from_records(d)
    else:
        df = None

    return df, mydate

def search_finder_files(mydate = None, user=None):
    '''
    Returns the files that are present in the disk at a given date.
    It also returns a message stating what date that was.

    TODO: This routine that looks in the disk, should look at the database and retrieve only the files which
    correspond to the privileges of the user.
    '''
    #If the date is specified, we will try to located the right file.
    #None will be returned if it does not exist.
    if ( not mydate is None):
        files = glob.glob(os.path.join("/scr2/sedm/phot", mydate, "finders/*ACQ*.jpg"))
        files = files + glob.glob(os.path.join("/scr2/sedm/phot", mydate, "finders/*ACQ*.png"))
        if len(files) == 0:
            files = None

    else:         
        curdate = datetime.datetime.utcnow()
        #Try to find the stat files up to 100 days before today's date.
        i = 0
        files = []
        while i < 100:
            newdate = curdate - datetime.timedelta(i)
            newdatedir = "%d%02d%02d"%(newdate.year, newdate.month, newdate.day)
            files = glob.glob(os.path.join("/scr2/sedm/phot", newdatedir, "finders/*ACQ*.jpg"))
            files = files + glob.glob(os.path.join("/scr2/sedm/phot", newdatedir, "finders/*ACQ*.png"))
            if len(files) > 0:
                mydate = newdatedir
                break

            i = i+1

    if not files is None:
        files.sort(reverse=True)
        d = {'filename':files}
        df = DataFrame.from_records(d)
    else:
        df = None

    return df, mydate

def search_phot_files_by_imtype(mydate = None, user=None):
    '''
    Returns the files that are present in the disk at a given date.
    It also returns a message stating what date that was.

    TODO: This routine that looks in the disk, should look at the database and retrieve only the files which
    correspond to the privileges of the user.
    '''
    #If the date is specified, we will try to located the right file.
    #None will be returned if it does not exist.

    filedic = {}

    files = glob.glob(os.path.join("/scr2/sedm/phot", mydate, "rc*.fits"))

    for f in files:
        imtype = fitsutils.get_par(f, "IMGTYPE").title()
        prev = filedic.get(imtype, [])
        path, fits = os.path.split(f)
        prev.extend([os.path.join(path, "pngraw", fits.replace(".fits", "_all.png") )])
        filedic[imtype] = prev

    return filedic


def get_requests_for_user(user_id, inidate, enddate):
    '''
    Obtains the DataFrame for the requests that were made:
        - By any member of the group where the user with user_id belongs to.
        - In the last 5 days
    '''
    request_query = ("""SELECT a.designator, o.name, o.ra, o.dec, r.inidate, r.enddate, r.priority, r.status, r.lastmodified, r.id 
                        FROM request r, object o, allocation a 
                        WHERE o.id = r.object_id AND a.id = r.allocation_id  
                            AND ( r.lastmodified >= DATE('%s') AND r.lastmodified <= DATE('%s') )
                            AND r.allocation_id IN
                           (SELECT a.id
                            FROM allocation a, groups g, usergroups ug, users u, program p
                            WHERE ug.user_id = u.id AND ug.group_id = g.id AND u.id = %d AND p.group_id = g.id AND a.program_id = p.id
                            ) ORDER BY r.lastmodified DESC;"""% (inidate, enddate, user_id))
    req = db.execute_sql(request_query)
    req = DataFrame(req, columns=['allocation', 'object', 'RA', 'DEC', 'start date', 'end date', 'priority','status', 'lastmodified', 'UPDATE'])
    if user_id in [189, 2, 20180523190352189]:
        req['UPDATE'] = req['UPDATE'].apply(convert_to_link)
        #req['object'] = req['object'].apply(convert_to_growth_link)
    else:
        req.drop(columns=['UPDATE', 'RA', 'DEC'])


    return req

def convert_to_link(reqid):

    return """<a href='view_request?id=%s'>+</a>""" % reqid

def convert_to_growth_link(reqid):

    return """<a href='http://skipper.caltech.edu:8080/cgi-bin/growth/view_source.cgi?name=%s'>%s</a>""" % (reqid, reqid)

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
    elif ("'" in username):
        username = username.split("'")[1]
        user_info = db.execute_sql("SELECT username, name, email, id FROM users WHERE username ='{0}'".format(username))
    else:
        userlower = username.lower()
        user_info = db.execute_sql("SELECT username, name, email, id FROM users WHERE LOWER(username) LIKE '%{0}%' OR LOWER(name) LIKE '%{1}%' OR LOWER(email) LIKE '%{2}%@%'".format(userlower, userlower, userlower))

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

    res = db.execute_sql(""" SELECT a.id, a.designator, p.designator, g.designator, a.time_allocated, a.time_spent
                            FROM allocation a, program p, groups g, usergroups ug
                            WHERE a.program_id = p.id AND p.group_id = g.id 
                            AND g.id = ug.group_id AND a.active is True AND ug.user_id = %d"""%(user_id))

    # create the dataframe and set the allocation names to be linked
    data = DataFrame(res, columns=['id', 'allocation', 'program', 'group', 'time allocated', 'time spent'])

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


def get_allocation_stats(user_id, inidate=None, enddate=None):
    """
    Obtains a list of allocations that belong to the user and 
    query the total allocated name and time spent for that allocation.

    If no user_id is provided, all active allocations are returned.
    """
    if (user_id is None):
        res = db.get_from_allocation(["designator", "time_allocated", "time_spent"], {"active":True})
        df = DataFrame(res, columns=["designator", "time_allocated", "time_spent"])

        alloc_hours = np.array([ta.total_seconds() / 3600. for ta in df["time_allocated"]])
        spent_hours = np.array([ts.total_seconds() / 3600. for ts in df["time_spent"]])
        free_hours = alloc_hours - spent_hours

        df = df.assign(alloc_hours=alloc_hours, spent_hours=spent_hours, free_hours=free_hours)

    else:
        if (inidate is None or enddate is None):
            res = db.execute_sql(""" SELECT a.designator, a.time_allocated, a.time_spent
                                    FROM allocation a, program p, groups g, usergroups ug
                                    WHERE a.program_id = p.id AND p.group_id = g.id 
                                    AND g.id = ug.group_id AND a.active is True AND ug.user_id = %d"""%(user_id))

            df = DataFrame(res, columns=["designator", "time_allocated", "time_spent"])

            alloc_hours = np.array([ta.total_seconds() / 3600. for ta in df["time_allocated"]])
            spent_hours = np.array([ts.total_seconds() / 3600. for ts in df["time_spent"]])
            free_hours = alloc_hours - spent_hours

            df = df.assign(alloc_hours=alloc_hours, spent_hours=spent_hours, free_hours=free_hours)


        else:
            res = db.execute_sql(""" SELECT DISTINCT a.id, a.designator, a.time_allocated
                                    FROM allocation a, program p, groups g, usergroups ug
                                    WHERE a.program_id = p.id AND p.group_id = g.id 
                                    AND g.id = ug.group_id AND a.active is True AND ug.user_id = %d;"""%(user_id))
            allocdes = []
            spent_hours = []
            alloc = []
            for ais in res:
                spent = db.get_allocation_spent_time(ais[0], inidate, enddate)
                allocdes.append(ais[1])
                spent_hours.append(int(spent)/3600.)
                alloc.append(ais[2])
            res = np.array([allocdes, alloc, spent_hours])

            df = DataFrame(res.T, columns=["designator", "time_allocated", "time_spent"])
            alloc_hours = np.array([ta.total_seconds() / 3600. for ta in df["time_allocated"]])
            free_hours = alloc_hours - spent_hours
            df = df.assign(alloc_hours=alloc_hours, spent_hours=spent_hours, free_hours=free_hours)


    df = df.sort_values(by=["alloc_hours"], ascending=False)

    alloc_names = df["designator"].values
    category = ["alloc_hours", "spent_hours", "free_hours"]


    data = {'allocations' : alloc_names}

    for cat in category:
        data[cat] = df.fillna(0)[cat]

    return data

if __name__ == "__main__":
    #get_request_by_id(20180416223941629)
    pass
