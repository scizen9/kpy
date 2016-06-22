# -*- coding: utf-8 -*-
"""
Created on Tue Mar  1 20:19:04 2016

@author: nadiablago
"""

import os
import shutil
import subprocess
import numpy as np
import pyfits as pf
import rcred
import datetime
import fitsutils

from matplotlib import pylab as plt


def run_sex(flist, mask=False, cosmics=True):
    
    d = os.path.dirname(flist[0])
    if d == "":
        d = "."
    os.chdir(d)

    #Create the directory where the sextracted images are going to go.
    sexdir = os.path.join(d, "sextractor")    
    if (not os.path.isdir(sexdir)):
        os.makedirs(sexdir)
        
    newlist = []
    for f in flist:
        try:
            f = os.path.abspath(f)
            if (mask):
                out = rcred.get_masked_image(f)
            else:
                out = f
                
            if (cosmics and (not fitsutils.has_par(out, "CRREJ") or fitsutils.get_par(out, "CRREJ") ==0)):
                out = rcred.clean_cosmic(out)

            cmd="sex -c %s/config/daofind.sex %s"%(os.environ["SEDMPH"], out) 
            subprocess.call(cmd, shell=True)
            print cmd
            newimage = os.path.join(sexdir, os.path.basename(f).replace(".fits", ".sex")) 
            shutil.move("image.sex", newimage)
            newlist.append(newimage)
        except IOError:
            print "IOError detected reading file",f
            pass
        
    return newlist
        
def analyse_sex(sexfileslist, plot=True, interactive=False):
    '''
    Analyses the sextractor filelist to determine the best focus.
   
	#   1 X_IMAGE                Object position along x                                    [pixel]
	#   2 Y_IMAGE                Object position along y                                    [pixel]
	#   3 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
	#   4 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
	#   5 MAG_BEST               Best of MAG_AUTO and MAG_ISOCOR                            [mag]
	#   6 MAGERR_BEST            RMS error for MAG_BEST                                     [mag]
	#   7 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
	#   8 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
	#   9 ELLIPTICITY            1 - B_IMAGE/A_IMAGE                                       
	#  10 BACKGROUND             Background at centroid position                            [count]
	#  11 FLAGS                  Extraction flags   
 
    returns: A tuple containing:
        1. - The best focus values as interpolated from the images.
        2. - The sigma whithin which we should look for a finer focus.
    '''
    sexfileslist = list(sexfileslist)    
    sexfileslist.sort()
 
    focpos = []
    fwhms = []
    std_fwhm = []
    for i, f in enumerate(sexfileslist):
        fits = f.replace("sextractor/", "").replace(".sex", ".fits")
        FF = pf.open(fits)
        pos= float(FF[0].header['focpos'])

        s = np.genfromtxt(f, comments="#")
        
        s = s[s[:,1]< 2000]

	#Only objects with FWHM less than 20 pixels...
        s = s[s[:,7] < 20]
        
        #Select round sources (ellipticity is 1-axis_ratio)
        s = s[s[:,8]<np.percentile(s[:,8], 30)]
        #Select bright magnitudes
        s = s[s[:,5]<np.percentile(s[:,5], 20)]
        print f, "number of sources", len(s)
 
        focpos.append(pos)
        fwhms.append(np.nanmean(s[:,7]*0.394))
	mad = np.median(np.abs(s[:,7] - np.nanmean(s[:,7])))/0.67448975019608171 * 0.394
        #std_fwhm.append(np.std(s[:,6]*0.394))
        std_fwhm.append(mad)
    
    focpos = np.array(focpos)
    fwhms = np.array(fwhms)
    std_fwhm = np.maximum(1e-5, np.array(std_fwhm))
    
    coefs = np.polyfit(focpos, fwhms, w=1/std_fwhm, deg=2)
    
    x = np.linspace(np.min(focpos), np.max(focpos), 100)
    p = np.poly1d(coefs)
    print "Best focus:%.2f"% x[np.argmin(p(x))], coefs,std_fwhm
    
    
    if (plot==True):
        plt.title("Best focus:%.2f"% x[np.argmin(p(x))])
        with open("/tmp/focus", "w") as f:
            f.write(str(focpos))
            f.write(str(fwhms))
        plt.errorbar(focpos, fwhms, yerr=std_fwhm, fmt="o")
        plt.plot(x, p(x), "-")
        plt.xlabel("Focus (mm)")
        plt.ylabel("FWHM (arcsec)")
        if (interactive):
            plt.show()
        else:
            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]),"focus_%s.png"%(datetime.datetime.utcnow()).strftime("%Y%m%d-%H:%M:%S")))
    return x[np.argmin(p(x))], coefs[0]
    

def analyse_image(sexfile):
    '''
    Analyses the sextractor filelist to determine the best focus.
    
    returns: A tuple containing:
        1. - Number of extracted sources.
        2. - FWHM in arcsecs.
	3. - Ellipticity.

	#   1 X_IMAGE                Object position along x                                    [pixel]
	#   2 Y_IMAGE                Object position along y                                    [pixel]
	#   3 ALPHA_J2000            Right ascension of barycenter (J2000)                      [deg]
	#   4 DELTA_J2000            Declination of barycenter (J2000)                          [deg]
	#   5 MAG_BEST               Best of MAG_AUTO and MAG_ISOCOR                            [mag]
	#   6 MAGERR_BEST            RMS error for MAG_BEST                                     [mag]
	#   7 FWHM_WORLD             FWHM assuming a gaussian core                              [deg]
	#   8 FWHM_IMAGE             FWHM assuming a gaussian core                              [pixel]
	#   9 ELLIPTICITY            1 - B_IMAGE/A_IMAGE                                       
	#  10 BACKGROUND             Background at centroid position                            [count]
	#  11 FLAGS                  Extraction flags   
 
    '''
    


    s = np.genfromtxt(sexfile, comments="#")

    if (s is None or s.ndim==0 or len(s)==0):
        return 0,0,0

    #Select sources inside of the cross
    x = s[:,0]
    y = s[:,1]
    s = s[((y<850)|(y>1125))*((x<885)|(x>1540))]

    # Select with good flags only.
    s = s[s[:,10]==0]

    nsources = len(s) 
    if (nsources > 0):
        print nsources
    else:
        return 0,0,0
    #Select round sources (ellipticity is 1-axis_ratio)
    s = s[s[:,8]<0.5]
    ellipticity = np.nanmedian(s[:,8])
    s = s[s[:,8]<0.25]

    #Select FWHM at least 3 pixels and lower than 6 arcsec
    s = s[ (s[:,7]>3)*(s[:,7]*0.394<6)]
    
    nsources = len(s) 
    if (nsources > 0):
        print nsources
    else:
        return 0,0,0
        
    #Select bright magnitudes
    s = s[s[:,4]<np.percentile(s[:,4], 20)]
       
    fwhm = np.nanmedian(s[:,7]*0.394)
    
    return nsources, fwhm, ellipticity
        
def get_focus(lfiles, plot=True):
    '''
    Receives a list of focus files and returns the best focus value.
    
    '''
    sexfiles = run_sex(lfiles)
    focus, sigma = analyse_sex(sexfiles, plot=plot)
    return focus, sigma
    
def get_image_pars(image):
    '''
    Returns a set of statistics for a given image.
    '''
    sexfiles = run_sex([image])
    pars = analyse_image(sexfiles[0])
    
    return pars
    
