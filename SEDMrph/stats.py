# -*- coding: utf-8 -*-
"""
Created on Sat May 21 10:26:37 2016

@author: nadiablago
"""

import glob, os
import recenter_ifu
import fitsutils
import coordinates_conversor as cc
import numpy as np
import sextractor 
import pyfits as pf
from matplotlib import pylab as plt

def compile_stats_pointing():
    ra = 0
    dec = 0
    
    out = open("/tmp/pointing", "w")
    out.write("#f, imtype, obj, jd, filter, radeg, decdeg, dra, ddec\n")
    myfiles = glob.glob("/scr2/sedm/phot/20160616/a_*[0-9].fits")
    myfiles.sort()
    for f in myfiles:
        #try:
        imtype = fitsutils.get_par(f, "IMGTYPE")
        newra = fitsutils.get_par(f, "OBJRA")
        newdec = fitsutils.get_par(f, "OBJDEC")
        newra, newdec = cc.hour2deg(newra, newdec)        
        myfilter =  fitsutils.get_par(f, "FILTER")
        if (imtype == "ACQUISITION" or imtype == "SCIENCE"):#: and np.round(ra, 2) != np.round(newra, 2) and np.round(dec, 2) != np.round(newdec, 2):
            obj  = fitsutils.get_par(f, "OBJECT")
            jd = fitsutils.get_par(f, "JD")
            status, dra, ddec = recenter_ifu.get_offset_center(f, plot=False, interactive=False)
        
            print f, ra, newra, imtype, jd, dra, ddec
            out.write("%s,%s,%s,%.2f,%s,%.5f,%.5f,%.2f,%.2f\n"%(f, imtype, obj, jd, myfilter, newra, newdec, dra, ddec))
        ra = newra
        dec = newdec
    
        #except:
        #    pass
    out.close()
    
def plot_stats_pointing():
    
    from matplotlib import pylab as plt
    
    d = np.genfromtxt("/tmp/pointing", dtype=None, delimiter=",", names=True)    
    
    
    plt.scatter(d["dra"], d["ddec"], c=d["jd"]-np.min(d["jd"]))
    cb = plt.colorbar()
    cb.set_label("JD- MIN(JD)")
    plt.xlabel("dRA [arcsec]")
    plt.ylabel("dDEC [arcsec]")
    plt.savefig("/tmp/pointing_errors.png")
    
    plt.xlim(-60,60)
    plt.ylim(-60,60)
    plt.savefig("/tmp/pointing_errors_60.png")
    plt.clf()
    
    d1 = d[d["imtype"]!="SCIENCE"]
    d = d[d["imtype"]=="SCIENCE"]
    dr = d[d["filter"]=="r"]
    dg = d[d["filter"]=="g"]
    di = d[d["filter"]=="i"]
    du = d[d["filter"]=="u"]
    
    filters = ["ACQ"]#, "r", "g", "i", "u" ]    
    
    dar = np.array([d1])#r, dg, di, du, d1])
    for i in np.arange(len(filters)):
        plt.figure(i+1)
        plt.scatter(dar[i]["dra"], dar[i]["ddec"], c=dar[i]["jd"]-np.min(d["jd"]))
        cb = plt.colorbar()
        plt.title(filters[i])
        cb.set_label("JD- MIN(JD)")
        plt.xlabel("dRA [arcsec]")
        plt.ylabel("dDEC [arcsec]")
        plt.savefig("/tmp/pointing_errors_%s.png"%filters[i])
        plt.clf()
    plt.close("all")


def get_sextractor_stats(files):
    
    files.sort()
    #sexfiles = sextractor.run_sex(files)
    sexfiles = [os.path.join(os.path.join(os.path.dirname(f), "sextractor"), os.path.basename(f).replace(".fits", ".sex")) for f in files]    
    
    sexfiles.sort()
    
    with open(os.path.join( os.path.dirname(files[0]), "stats/stats.log"), "w") as out:
        for i, f in enumerate(files):
            print sexfiles[i]
            if not os.path.isfile(sexfiles[i]):
                sf =  sextractor.run_sex([f])
            else:
                sf = sexfiles[i]
            print f
            hd = pf.open(files[i])[0].header
            try:
                jd = hd["JD"]
                ns, fwhm, ellipticity = sextractor.analyse_image(sf)
                out.write("%s,%.3f,%d,%.2f,%.3f\n"%(os.path.abspath(f),jd,ns,fwhm,ellipticity))
            except:
                pass
            
def plot_stats(statfile):
    s = np.genfromtxt(statfile, delimiter=",", dtype=None)
    s.sort(order="f1")
    s = s[s["f2"]>0]
    f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2, sharex='col')
    ax1.plot(s["f1"], s["f2"], ".-")
    ax1.set_title('Number of bright sources extracted')
    ax2.plot(s["f1"], s["f3"], ".-")
    ax2.set_title('FWHM [arcsec]')
    ax3.plot(s["f1"], s["f4"], ".-")
    ax3.set_title('Ellipticity')
    #ax4.plot(x, 2 * y ** 2 - 1, color='r')
    plt.show()