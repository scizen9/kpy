# -*- coding: utf-8 -*-
"""
Created on Thu Jun 15 10:50:10 2017

@author: nadiablago
"""

from __future__ import print_function 
import os, glob, shutil
import numpy as np
import argparse
from SEDMrph import fitsutils
from astropy.time import Time
import datetime
from astropy.io import fits
import archive_pop_sedmdb

#cur_dir = os.getcwd()
#os.chdir('/home/sedm/kpy/SEDMDb')
sdb = archive_pop_sedmdb.SedmDB(dbname='sedmdbtest', host='localhost')
#os.chdir(cur_dir)

def ra_to_decimal(ra):
    hms = ra.split(':')
    return float(hms[0])*15+float(hms[1])/4+float(hms[2])/240


def dec_to_decimal(dec):
    dms = dec.split(':')
    return float(dms[0])+float(dms[1])/60+float(dms[2])/3600

def fill_par_dic_obs(fitsfile):
    '''
    Parses the fits file to fill the needed parameters for the database
    '''
    pardic_obs = {
    'id': 0,
    'object_id': 0, 
    'request_id': 0, 
    'mjd': 0.0, 
    'airmass': 0.0,
    'exptime': 0.0, 
    'fitsfile': '', 
    'lst': '', 
    'ra': 0, 
    'dec': 0, 
    'tel_ra': '',
    'tel_dec': '', 
    'tel_az': 0, 
    'tel_el': 0, 
    'tel_pa': 0, 
    'ra_off': 0,
    'dec_off': 0, 
    'imtype': '', 
    'camera': ''
            }
            
    f = fits.open(fitsfile)
    
    for k in pardic_obs.keys():
        if k in f[0].header:
            pardic_obs[k] = f[0].header[k]
    pardic_obs['mjd'] = f[0].header['JD'] - 2400000.5
    pardic_obs['ra'] = ra_to_decimal(f[0].header['RA'])
    pardic_obs['dec'] = dec_to_decimal(f[0].header['DEC'])
    pardic_obs['camera'] = f[0].header['CAM_NAME']   
    
        
    #Rtrieve the old obs_id from the header
    obs_id = f[0].header['obs_id']
    pardic_obs['id'] = obs_id

    return pardic_obs


def create_user(f):
    '''
    Reads the email encoded in the header and assigns it to a user.
    It tries to locate a user in te DB with that email.
    If not, a new user is created.
    '''
    email = fitsutils.get_par(f, "EMAIL")
    res = sdb.get_from_users(["id"], {"email":email})
    if len(res) == 1:
        return res[0]
    else:
        username = email.split("@")[0]
        sdb.add_user({'username':username,
                      'name':'Unknown',
                      'email':email,
                      'password':username})
        res = sdb.get_from_users(["id"], {"email":email})
        
        sdb.add_group({'designator': ''})
        return res[0]
    
def create_program(f):
    '''
    Creates a program registry in the program table.
    
    '''
    program_id = fitsutils.get_par(f, "P60PRID").uppercase()
    program_name = fitsutils.get_par(f, "P60PRNM")
    program_pi = fitsutils.get_par(f, "P60PRPI")
    target_id = fitsutils.get_par(f, "TARGID")
    
    pgrogid = sdb.get_from_program(["id"], {"descriptor":program_id})
    
    if (len(res)==0):
        print ("Program with program id %s not found"%program_id)
    else:
        #TODO - create a new program with that ID
        #progid
        pass
    
    
    
                            
def create_cal_request(files, inidate, enddate, caltype):
    '''
    Creates a default calibration request assigned to sedmcal, the calibration user.
    It also looks at each file and creates an atomic request associated to each observation.
    '''
    caltype_obj = {'bias': 3, 'twilight':4, 'dome':5, 'focus':6, 'test':7, 'test2':8, 'test3':9}
    obj_id = caltype_obj[caltype]
    
    #First check that there is no request for that night
    res = sdb.get_from_request(['object_id'],
    {'object_id': obj_id, 'user_id':32, 'allocation_id':1, 'inidate':inidate, 'enddate':enddate},
    {'inidate':'>=', 'enddate':'<='})
    
    #If there is no request, we add it.
    if len(res) == 0:
        exptime = fitsutils.get_par(files[0], "EXPTIME")
        pardic = {'object_id':obj_id,
                  'user_id':32,
                  'allocation_id':1,
                  'exptime':'{0, %d}'%exptime,
                  'priority':.1,
                  'inidate': inidate, #('year-month-day') (start of observing window),
                  'enddate': enddate, #('year-month-day') (end of observing window),
                  'nexposures':'{0, %s, 0, 0, 0}' % len(files)}
        req = sdb.add_request(pardic)
        if req >= 0:
            return req[0]
        else:
            print(req[1])
        
def create_obs(reqid, files, inidate, enddate, caltype):
    '''
    For each file, creates the observation associated with the file.
    
    '''
    caltype_obj = {'bias': 3, 'twilight':4, 'dome':5, 'focus':6, 'test':7, 'test2':8, 'test3':9}
    obj_id = caltype_obj[caltype]

    for f in files:
                
        jd_init = fitsutils.get_par(f, "JD")
        jd_end = jd_init + fitsutils.get_par(f, "EXPTIME")/(3600*24.)

        inidate = Time(jd_init, format='jd').iso
        enddate = Time(jd_end, format='jd').iso
                        
        pardic_obs = fill_par_dic_obs(f)
        pardic_obs["object_id"] = obj_id
        pardic_obs["request_id"] = reqid
        pardic_obs['fitsfile'] = f
        pardic_obs['imtype'] = caltype
        
        #Also add the observation

        obs = sdb.add_observation(pardic_obs)
        if obs[0] < 0:
            print(obs[1])


    
def log_calibrations(lfiles, caltype="test"):
    '''
    Logs the calibrations.
    
    These are the default objects of calibration.
    
    (id, marshal_id, name, ra, dec) (3, 0, 'bias', 0,0);
    (id, marshal_id, name, ra, dec) values (4, 1, 'twilight', 0,0);
    (id, marshal_id, name, ra, dec) values (5, 2, 'dome', 0,0);
    (id, marshal_id, name, ra, dec) values (6, 3, 'focus', 0,0);
    (id, marshal_id, name, ra, dec) values (7, 4, 'test', 0,0);
    (id, marshal_id, name, ra, dec) values (8, NULL, 'test2', 0,0);
    (id, marshal_id, name, ra, dec) values (9, NULL, 'test3', 0,0);

    '''
    
    #Retrieve the date of the first file
    lfiles.sort()
    
    jd_init = fitsutils.get_par(lfiles[0], "JD")
    jd_end = fitsutils.get_par(lfiles[-1], "JD")
    
    inidate = Time(jd_init, format='jd').iso
    enddate = Time(jd_end, format='jd').iso
    
    #First make sure that we have a Calibration Request to which assign all the registers.
    reqid = create_cal_request(lfiles, inidate, enddate, caltype)
    
    #Then associate all the files to that request id.
    create_obs(reqid, lfiles, inidate, enddate, caltype)

    sdb.update_request({'id': reqid, 'STATUS': 'COMPLETED'})
    
    
    
def log_science(lfiles, scitype):
    '''
    Logs the science files.
    '''
    
    for f in lfiles:
        #Initial time of the request is the starting JD.
        jd_init = fitsutils.get_par(f, "JD")
        #The end time of the request is the starting point plus the exposure plus 60s overhead.
        jd_end = fitsutils.get_par(f, "JD") + (fitsutils.get_par(f, "EXPTIME") + 60)/(24*3600.) 
        
        inidate = Time(jd_init, format='jd').iso
        enddate = Time(jd_end, format='jd').iso
    
        #First make sure that we have an Observation Request to which assign all the registers.
        #For observations with the same object name we will consider that they are associated with the same request id.
        reqid = create_cal_request(inidate, enddate, caltype)
    
    if scitype == "GUIDER":
        #log guider
        pass
    elif scitype== "SCIENCE":
        #log science file
        pass
    elif scitype== "ACQUISITON":
        #log acquisition file
        pass
    #Call fill_part_dic_obs 
    #to retrieve all the dicitonary fields needed to fill in the OBSERVATION table
    
    #Add the observation to the DB with the right user/group/request
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

        Parses a directory with one night photometric data.
        Obtains the relevant values from the headers and populates the DB.
        
        %run populate_db_from_archive.py -d PHOTDIR 
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-d', '--photdir', type=str, help='Directory containing the science fits for the night.', default=None)

    args = parser.parse_args()
    
    photdir = args.photdir
    

    myfiles = {
        "ACQUISITION":[],
        "BIAS":[],
        "DOME":[],
        "FOCUS":[],
        "GUIDER":[],
        "SCIENCE":[],
        "TWILIGHT":[]}
        
    mydir = os.path.abspath(photdir)
    #Gather all RC fits files in the folder with the keyword IMGTYPE=SCIENCE
    for f in glob.glob(os.path.join(mydir, "rc*fits")):
        try:
            if (fitsutils.has_par(f, "IMGTYPE")):
                imgtype = fitsutils.get_par(f, "IMGTYPE")
                myfiles[imgtype].append(f)
        except:
            print("problems opening file %s" %f)
         
    log_calibrations(myfiles["BIAS"], caltype="bias")
    #log_calibrations(myfiles["DOME"], caltype="dome")
    #log_calibrations(myfiles["FOCUS"], caltype="focus")
    #log_calibrations(myfiles["TWILIGHT"], caltype="twilight")
    
    #log_science(myfiles["SCIENCE"], "science")
    #log_science(myfiles["GUIDER"], "guider")
    #log_science(myfiles["ACQUISITION"], "acquisition")
