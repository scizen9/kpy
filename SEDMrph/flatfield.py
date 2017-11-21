# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 15:34:47 2016

@author: nadiablago
"""
import glob
import fitsutils
import rcred
import pyfits as pf
import numpy as np
import datetime
import time_utils
import coordinates_conversor as cc
from matplotlib import pylab as plt
from scipy.optimize import curve_fit

        
def get_flats_counts(directory):

    flatlist = []
    corners = {
    "g" : [1, 910, 1, 900],
    "i" : [1, 910, 1060, 2045],
    "r" : [1040, 2045, 1015, 2045],
    "u" : [1030, 2045, 1, 900]
    }
    
    for f in glob.glob(directory + "/rc*fits"):
        if fitsutils.has_par(f, "IMGTYPE") and fitsutils.get_par(f, "IMGTYPE") == "TWILIGHT":
            flatlist.append(f)
            
    if (len(flatlist)==0):
        print "No suitable twilight flats found in directory: %s"%directory
        return
        
    counts = {"u":[], "g":[], "r":[], "i":[]}
    sun_decs = {"u":[], "g":[], "r":[], "i":[]}
    
    for f in flatlist:
        print f
        print fitsutils.get_par(f, "JD")
        
        bias = rcred.get_overscan_bias_rc(f)
        exptime = fitsutils.get_par(f, "EXPTIME")
        sunsettime = fitsutils.get_par(f, "SUNSET")
        utc = time_utils.jd2utc(fitsutils.get_par(f, "JD"))
        st = datetime.datetime.strptime(sunsettime, "%H:%M")

        
        elapsed = 3600*(utc.hour - st.hour) + 60*(utc.minute - st.minute) + (utc.second - st.second)
        
        if (elapsed > 3000):
            continue
        print elapsed, utc, st
        
        
        data = pf.open(f)[0].data
        for band in corners.keys():
            c = corners[band]
            sf = data[c[0]:c[1], c[2]:c[3]]
            if (np.percentile(sf, 90) < 55000):
                counts[band].append( (np.percentile(sf, 90)-bias)/exptime)
                sun_decs[band].append(elapsed)        
    
    


    for band in corners.keys():
        coefs = np.polyfit(sun_decs[band], np.log10(counts[band]), deg=1, w=1./np.sqrt(np.log10(counts[band])))
        p = np.poly1d(coefs)
        x = np.linspace(np.min(sun_decs[band]), np.max(sun_decs[band]), 1000)
        plt.plot(sun_decs[band], np.log10(counts[band]), "o", label=band)
        plt.plot(x, p(x), label="Model "+band)
    plt.xlabel("Elapsed second since Sunset")
    plt.ylabel("Counts/s")
    plt.legend()
    plt.show()
            
