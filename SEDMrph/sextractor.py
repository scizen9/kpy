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

from matplotlib import pylab as plt


def run_sex(flist):
    
    d = os.path.dirname(flist[0])
    os.chdir(d)

    #Create the directory where the sextracted images are going to go.
    sexdir = os.path.join(d, "sextractor")    
    if (not os.path.isdir(sexdir)):
        os.makedirs(sexdir)
        
    newlist = []
    for f in flist:
        f = os.path.abspath(f)
        masked = rcred.get_masked_image(f)
        cmd="sex -c %s/config/daofind.sex %s"%(os.environ["SEDMPH"], masked) 
        subprocess.call(cmd, shell=True)
        print cmd
        newimage = os.path.join(sexdir, os.path.basename(f).replace(".fits", ".sex")) 
        shutil.move("image.sex", newimage)
        newlist.append(newimage)
        
    return newlist
        
def analyse_sex(sexfileslist, plot=True, interactive=False):
    '''
    Analyses the sextractor filelist to determine the best focus.
    
    returns: A tuple containing:
        1. - The best focus values as interpolated from the images.
        2. - The sigma whithin which we should look for a finer focus.
    '''
    
    focpos = []
    fwhms = []
    std_fwhm = []
    for i, f in enumerate(sexfileslist):
        fits = f.replace("sextractor/", "").replace(".sex", ".fits")
	print fits
        FF = pf.open(fits)
        pos= float(FF[0].header['focpos'])

        s = np.genfromtxt(f, comments="#")
        
        s = s[s[:,1]< 2000]
        
        #Select round sources (ellipticity is 1-axis_ratio)
        s = s[s[:,7]<np.percentile(s[:,7], 15)]
        #Select bright magnitudes
        s = s[s[:,2]<np.percentile(s[:,2], 15)]
        print "number of sources", len(s)
 
        focpos.append(pos)
        fwhms.append(np.nanmean(s[:,6]*0.394))
        std_fwhm.append(np.std(s[:,6]*0.394))
    
    focpos = np.array(focpos)
    fwhms = np.array(fwhms)
    std_fwhm = np.maximum(1e-5, np.array(std_fwhm))
    
    coefs = np.polyfit(focpos, fwhms, w=1/std_fwhm, deg=2)
    
    x = np.linspace(np.min(focpos), np.max(focpos), 100)
    p = np.poly1d(coefs)
    print "Best focus:%.2f"% x[np.argmin(p(x))], coefs,std_fwhm, p(x)
    
    
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
            plt.savefig(os.path.join(os.path.dirname(sexfileslist[0]),"focus.png"))
    return x[np.argmin(p(x))], coefs[0]
    

def analyse_image(sexfile):
    '''
    Analyses the sextractor filelist to determine the best focus.
    
    returns: A tuple containing:
        1. - The best focus values as interpolated from the images.
        2. - The sigma whithin which we should look for a finer focus.
    '''
    


    s = np.genfromtxt(sexfile, comments="#")
    #Select round sources (ellipticity is 1-axis_ratio)
    s = s[s[:,7]<0.1]

    #Select FWHM at least 1.2 arcsec and lower than 6
    s = s[ (s[:,6]*0.394>1.2)*(s[:,6]*0.394<6)]
    
    nsources = len(s)
    
    if (nsources > 0):
        print nsources
    else:
        return 0,0,0
        
    #Select bright magnitudes
    s = s[s[:,2]<np.percentile(s[:,2], 20)]
       
    fwhm = np.nanmedian(s[:,6]*0.394)
    ellipticity = np.nanmedian(s[:,7]*0.394)
    
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
    
