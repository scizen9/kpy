# -*- coding: utf-8 -*-
"""
Created on Wed Apr  6 16:13:31 2016

@author: nblago
"""

from numpy import *
from math import *
from types import *
import time_utils

def update_fremling():
    racen=256.9935
    deccen=34.3311
    issub='t'
    refsys='SDSS'
    filter='g'
    exptime=240
    utdate='2016-04-06T 12:41:45.6000'
    mag=99
    magerr=99
    limmag=18.7
    db = pg.DB(dbname='ptftransient',user='iptf_user',passwd='followup',host='yupana.caltech.edu')
    getsrcidquery = "SELECT id from sources where name='16yf'"
    result = db.query(getsrcidquery)
    for row in result.dictresult():
             srcid =  int(row['id'])
    query3 = "SELECT id from phot WHERE sourceid=%d and instrumentid=64 and obsdate='%s' and filter='%s' and reducedby='Fremling Automated Pipeline SEDM';" % (srcid,utdate,filter)
    result=db.query(query3)
    print query3
    if len(result.dictresult()) == 0:
             query2 = "INSERT INTO phot (sourceid,programid,instrumentid,ra,dec,obsdate,exptime,filter,mag,emag,limmag,issub,refsys,observer,reducedby) VALUES "
             query2 += "(%d,1,64,%f,%f,'%s',%f,'%s',%f,%f,%f,'%s','%s','SEDMachine','Fremling Automated Pipeline SEDM');" % (srcid,racen,deccen,utdate,exptime,filter,mag,magerr,limmag,issub,refsys)
             print query2
             db.query(query2)
    else:
        for row in result.dictresult():
            photid =  int(row['id'])
            query2 = "DELETE from phot WHERE id=%d;" % (photid)
            print query2
            db.query(query2)
            
            query2 = "INSERT INTO phot (sourceid,programid,instrumentid,ra,dec,obsdate,exptime,filter,mag,emag,limmag,issub,refsys,observer,reducedby) VALUES "
            query2 += "(%d,1,64,%f,%f,'%s',%f,'%s',%f,%f,%f,'%s','%s','SEDMachine','Fremling Automated Pipeline SEDM');" % (srcid,racen,deccen,utdate,exptime,filter,mag,magerr,limmag,issub,refsys)
            print query2
            db.query(query2)
         
         
def update_phot_blagorodnova(image):
    '''
    Updates the DB with the aperture photometry extrated from the fits and the interpolated zeropoint.
    '''
    racen=fitsutils.get_par(image, "OBJRA")
    deccen=fitsutils.get_par(image, "OBJRA")
    issub='f'
    refsys='SDSSinterpolated'
    filter=fitsutils.get_par(image, "FILTER")
    exptime=fitsutils.get_par(image, "EXPTIME")
    utdate=time_utils.jd2utc(fitsutils.get_par(image, "JD"), string=True)
    mag = np.round(fitsutils.get_par(image, "APPMAG"), 3)
    magerr = np.round(fitsutils.get_par(image, "APPMAGER"), 3)
    limmag = np.round(fitsutils.get_par(image, "ZEROPT"), 3)
    name =   fitsutils.get_par(image, "NAME").replace("PTF", "")
    observer = 'SEDMachine'
    reducedby = 'Blagorodnova Automated Pipeline SEDM'
    
    db = pg.DB(dbname='ptftransient',user='iptf_user',passwd='followup',host='yupana.caltech.edu')
    getsrcidquery = "SELECT id from sources where name='%s'"%name
    result = db.query(getsrcidquery)
    for row in result.dictresult():
             srcid =  int(row['id'])
    query3 = "SELECT id from phot WHERE sourceid=%d and instrumentid=64 and obsdate='%s' and filter='%s' and reducedby='Blagorodnova Automated Pipeline SEDM';" % (srcid,utdate,filter)
    result=db.query(query3)
    print query3
    
    if len(result.dictresult()) == 0:
             query2 = "INSERT INTO phot (sourceid,programid,instrumentid,ra,dec,obsdate,exptime,filter,mag,emag,limmag,issub,refsys,observer,reducedby) VALUES "
             query2 += "(%d,1,64,%f,%f,'%s',%f,'%s',%f,%f,%f,'%s','%s',%s,%s);" % (srcid,racen,deccen,utdate,exptime,filter,mag,magerr,limmag,issub,refsys,observer,reducedby)
             print query2
             db.query(query2)
    else:
        for row in result.dictresult():
            photid =  int(row['id'])
            query2 = "DELETE from phot WHERE id=%d;" % (photid)
            print query2
            db.query(query2)
            
            query2 = "INSERT INTO phot (sourceid,programid,instrumentid,ra,dec,obsdate,exptime,filter,mag,emag,limmag,issub,refsys,observer,reducedby) VALUES "
            query2 += "(%d,1,64,%f,%f,'%s',%f,'%s',%f,%f,%f,'%s','%s',%s,%s);" % (srcid,racen,deccen,utdate,exptime,filter,mag,magerr,limmag,issub,refsys,observer,reducedby)
            print query2
            db.query(query2)