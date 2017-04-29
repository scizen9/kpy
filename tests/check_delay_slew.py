# -*- coding: utf-8 -*-
"""
Created on Fri Mar 24 12:57:27 2017

@author: nblago
"""

import numpy as np
from matplotlib import pylab as plt
from astropy.coordinates import SkyCoord
import astropy.units as u

def get_lag(f):
    '''
    Input file obtained by running:
    gethead JD EXPTIME OBJECT EXPTIME RA DEC IMGTYPE 2017*/rc*fits | grep SCIENCE > ~/Documents/sedm_stats_slew
    
    on pharos
    
    '''
    a = np.genfromtxt(f, dtype=None)
    initime = a["f1"][:-1]
    endtime = a["f1"][1:]
    exptime = a["f2"][:-1]
    
    ra1 = a["f7"][:-1]
    dec1 = a["f8"][:-1]
    ra2 = a["f7"][1:]
    dec2 = a["f8"][1:]
    
    print ra1
    
    c1 = SkyCoord(ra1, dec1, unit=(u.hourangle, u.deg), obstime="J2000")
    c2 = SkyCoord(ra2, dec2, unit=(u.hourangle, u.deg), obstime="J2000")    
    
    angle = c1.separation(c2)
    
    lag = (endtime - (initime + (exptime+40)/(24.*3600))) * 24*3600
    same_obj = a["f3"][:-1] == a["f3"][1:]
    
    same_obj = same_obj[lag<300]
    lag = lag[lag<300]
    angle = angle[lag<300]
    
    plt.hist(lag[same_obj], range=(0,200), bins=50, label="Same object", normed=True, alpha=0.8)
    plt.hist(lag[~same_obj], range=(0,200), bins=50, label="New object", normed=True, alpha=0.8)
    
    plt.title("Median same object %.1fs \n Different object %.1fs"%(np.median(lag[same_obj]), np.median(lag[~same_obj])))
    plt.xlabel("Elapsed seconds after readout")
    plt.legend()
    
    plt.figure()
    plt.scatter(angle.value[~same_obj], lag[~same_obj]+40)
    plt.xlabel("Angle between slews [deg]")
    plt.ylabel("Total elapsed time [s]")
    plt.show()
    
    plt.figure()
    plt.scatter(angle.value[(angle.value<13/60.) * ~same_obj]*60, lag[~same_obj* (angle.value<13/60.)] + 40)
    plt.xlabel("Angle between slews [arcmin]")
    plt.ylabel("Total elapsed time [s]")
    plt.show()