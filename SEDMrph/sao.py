# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 04:01:07 2016

@author: nadiablago
"""

import ephem
import datetime
import coordinates_conversor
import os, subprocess
import numpy as np

def get_sao(radius=2000):
    
    '''
    Uses the current time and the latitude of Palomar to find the best SAO stars at zenith.
    
    Palomar.lon, Palomar.lat = '243.1361', '33.3558'    
    
    '''
    
    d = datetime.datetime.now()
    utc = d.utcnow()
    
    #Get reasonably high target
    ra = 15*((utc.hour+0.01)%24) + 15*(utc.minute/60.)
    dec = 40
    
    hra, hdec  =  coordinates_conversor.deg2hour(ra, dec)

    sao = get_sao_rec(hra, hdec, radius)
    while (len(sao) == 0):
        radius = radius+1000
        sao = get_sao_rec(hra, hdec, radius)
    print "Found %d"%len(sao)
    if np.ndim(sao) > 1:
        np.random.shuffle(sao)
    
    print "Returning ", sao[0][1], sao[0][2]
    
    saora, saodec  =  coordinates_conversor.hour2deg(sao[0][1], sao[0][2])
    
    return "SAO%s"%(sao[0][0]), saora, saodec 

def get_sao_rec(hra, hdec, radius):
    '''
    Runs the scat program to get an SAO star.
    '''
    print "Looking for an SAO at %s %s with radius %d"%(hra, hdec,radius)
    sao_path = os.environ["SAOCAT"]
    os.chdir(sao_path)
    cmd = "scat -c sao  %s %s J2000 -r %d -mx 7,9 > /tmp/sao"%(hra, hdec, radius)
    subprocess.call(cmd, shell=True)
    
    try:
        sao = np.genfromtxt("/tmp/sao", dtype="str")
        sao = np.array(sao, ndmin=1)
        print sao
    except IOError:
        sao = np.array([])
    
    return sao
#return saostars[]
#p = Popen(cmd.split(), stdout=PIPE, stdin=PIPE, stderr=STDOUT)
#es = p.communicate()

