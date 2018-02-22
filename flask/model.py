import numpy as np
import json
from datetime import datetime, timedelta
import os, glob
from SEDMDb.SedmDb import SedmDB 
from pandas import DataFrame

db = SedmDB(host='localhost', dbname='sedmdb')

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
        
        files = glob.glob(s) + glob.glob(s2) + glob.glob(s3)

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
            files = glob.glob(s) + glob.glob(s2) + glob.glob(s3)
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
