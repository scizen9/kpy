# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 13:11:52 2015

@author: nadiablago
"""
import fitsutils
import os, glob, shutil
import numpy as np
from pyraf import iraf 
import pyfits as pf
from matplotlib import pylab as plt
import subprocess

def create_masterbias(biasdir=None, channel='rc'):
    '''
    Combines slow and fast readout mode biases for the specified channel.
    '''
    
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    
    if (biasdir == None) or biasdir=="": biasdir = "."
        
    outs = "Bias_%s_slow.fits"%channel
    outf = "Bias_%s_fast.fits"%channel

    if (os.path.isfile(os.path.join(biasdir,outs)) and os.path.isfile(os.path.join(biasdir,outf))):
        print "Master Bias exists!"
        return
    else:
        print "Starting the Master Bias creation!"

    os.chdir(biasdir)        
        
    lfastbias = []
    lslowbias = []
    
    #Select all filts that are Bias with same instrument
    for f in glob.glob("*fits"):
        try:
            if (channel == fitsutils.get_par(f, "CHANNEL") and "BIAS" in str.upper(fitsutils.get_par(f, "OBJECT")) ):
                if (fitsutils.get_par(f, "ADCSPEED")==2):
                    lfastbias.append(f)
                else:
                    lslowbias.append(f)
        except:
            pass
                
    print "Files for bias SLOW mode: ", lslowbias
    print "Files for bias FAST mode: ", lfastbias
    
    if len(lfastbias) > 0:
        bfile_fast ="lbias_fast_"+channel
        np.savetxt(bfile_fast, np.array(lfastbias), fmt="%s")
        if (os.path.isfile("Bias_stats_fast")): os.remove("Bias_stats_fast")
        iraf.imstat("@"+bfile_fast, Stdout="Bias_stats_fast")
        
        st = np.genfromtxt("Bias_stats_fast", names=True, dtype=None)
        print st
        
        iraf.imcombine(input = "@"+bfile_fast, \
                    output = outf, \
                    combine = "median",\
                    scale = "mode")
        os.remove(bfile_fast)

    if len(lslowbias) > 0:

        bfile_slow ="lbias_slow_"+channel
        np.savetxt(bfile_slow, np.array(lslowbias), fmt="%s")
        if (os.path.isfile("Bias_stats_slow")): os.remove("Bias_stats_slow")
        iraf.imstat("@"+bfile_slow, Stdout="Bias_stats_slow")
        
        st = np.genfromtxt("Bias_stats_slow", names=True, dtype=None)
        print st
        
        iraf.imcombine(input = "@"+bfile_slow, \
                    output = outs, \
                    combine = "median",\
                    scale = "mode")
        os.remove(bfile_slow)
        
def solve_astrometry(img, radius=6.5, with_pix=True):
    '''
    img: fits image where astrometry should be solved.
    radius: radius of uncertainty on astrometric position in image.
    '''

    ra = fitsutils.get_par(img, 'RA')
    dec = fitsutils.get_par(img, 'DEC')
    print "Solving astrometry on field with (ra,dec)=", ra, dec
    
    astro = "a_" + img
    cmd = "solve-field --ra %s --dec %s --radius %.4f -p --new-fits %s \
      -W none -B none -P none -M none -R none -S none --overwrite %s"%(ra, dec, radius, astro, img)
    if (with_pix):
        cmd = cmd + "--scale-units arcsecperpix  --scale-low 0.375 --scale-high 0.425"
    print cmd

    subprocess.call(cmd, shell=True)
    
    #Cleaning after astrometry.net
    if (os.path.isfile(img.replace(".fits", ".axy"))):
        os.remove(img.replace(".fits", ".axy"))
    if (os.path.isfile(img.replace(".fits", "-indx.xyls"))):
        os.remove(img.replace(".fits", "-indx.xyls"))
    if (os.path.isfile("none")):
        os.remove("none")
        
    return astro
    
    
def slice_rc(img):
    '''
    Slices the Rainbow Camera into 4 different images and adds the 'filter' keyword in the fits file.
    '''
    fname = os.path.basename(img)
    fdir = os.path.dirname(img)    
    
    # Running IRAF
    iraf.noao(_doprint=0)
    
    '''f = pf.open(img)
    h = f[0].header
    i = f[0].data[1000:-1,0:920]
    g = f[0].data[0:910,0:900]
    r = f[0].data[1035:2046, 1050:2046]
    u = f[0].data[0:910,1050:2046]'''
    
    
    corners = {
    "i" : [1, 910, 1, 900],
    "g" : [1, 910, 1060, 2045],
    "r" : [1040, 2045, 1015, 2045],
    "u" : [1030, 2045, 1, 900]
    }
    
    
    filenames = []
        
    for i, b in enumerate(corners.keys()):
        name = fname.replace(".fits", "_%s.fits"%b)
        iraf.imcopy("%s[%d:%d,%d:%d]"%(img, corners[b][0], corners[b][1], corners[b][2], corners[b][3]), name)
         
        fitsutils.update_par(name, 'filter', b)

        filenames.append(name)
    
    return filenames

def clean_cosmic(f, name):
    '''
    From lacosmic.
    '''
    import cosmics
    
    
    g = fitsutils.get_par(f, "GAIN")
    rn = fitsutils.get_par(f, "RDNOISE")
    array, header = cosmics.fromfits(f)
    
    try:
        c = cosmics.cosmicsimage(array, gain=g, readnoise=rn, sigclip = 8.0, sigfrac = 0.3, satlevel = 64000.0)
        c.run(maxiter = 10)
        out = f.replace('.fits',  '_clean.fits')
    
        cosmics.tofits(out, c.cleanarray, header)
        
        os.remove(f)
    except:
        pass
    

    
                
def get_bias_rc(img):
    f = pf.open(img)
    bias = np.median(f[0].data[900:1100,900:1100].flatten())
    
    return bias
    
    
def create_masterflat(flatdir=None, biasdir=None, channel='rc'):
    '''
    Creates a masterflat from both dome flats and sky flats if the number of counts in the given filter
    is not saturated and not too low (between 1500 and 40000). 
    '''
    
    
    if (flatdir == None or flatdir==""): flatdir = "."
        
    if (biasdir == None or biasdir==""): biasdir = "."
        
    os.chdir(flatdir)
    
    if (len(glob.glob("Flat_%s*norm.fits"%channel)) == 4):
        print "Master Flat exists!"
        return 
    else:
        print "Starting the Master Flat creation!"

    bias_slow = "Bias_%s_fast.fits"%channel
    bias_fast = "Bias_%s_fast.fits"%channel
    
    if (not os.path.isfile(bias_slow) and not os.path.isfile(bias_fast) ):
        create_masterbias(biasdir)
     
    lsflat = []
    lfflat = []
    
    #Select all filts that are Flats with same instrument
    for f in glob.glob("*fits"):
        #try:
        if fitsutils.has_par(f, "OBJECT"):
            obj = str.upper(fitsutils.get_par(f, "OBJECT"))
        else:
            continue
        
        if ( ("DOME" in  obj or "FLAT" in obj) and (channel == fitsutils.get_par(f, "CHANNEL"))):
            if (fitsutils.get_par(f, "ADCSPEED")==2):
                lfflat.append(f)
            else:
                lsflat.append(f)
        #except:
        #    print "Error with retrieving parameters for file", f
        #    pass
                
    print "Files for slow flat", lsflat
    print "Files for fast flat", lfflat
    
    fsfile ="lflat_slow_"+channel
    np.savetxt(fsfile, np.array(lsflat), fmt="%s")
    fffile ="lflat_fast_"+channel
    np.savetxt(fffile, np.array(lfflat), fmt="%s")



    # Running IRAF
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    
    #Remove bias from the flat
    if len(lsflat) >0:
        iraf.imarith("@"+fsfile, "-", bias_slow, "b_@"+fsfile)
    
    if len(lfflat) >0:
        iraf.imarith("@"+fffile, "-", bias_fast, "b_@"+fffile)    
    
    #Slices the flats.
    debiased_flats = glob.glob("b_*.fits")
    for f in debiased_flats:
        print "Slicing file", f
        slice_rc(f)
        #Remove the un-sliced file
        os.remove(f)
        
    #Selects the ones that are suitable given the number of counts and combines them.
    bands = ['u', 'g', 'r', 'i']
    for b in bands:
        out = "Flat_%s_%s.fits"%(channel, b)
        out_norm = out.replace(".fits","_norm.fits")

        if (os.path.isfile(out_norm)):
            print "Master Flat for filter %s exists. Skipping..."%b
            continue
        
        lfiles = []
        for f in glob.glob('b_*_%s.fits'%b):
            d = pf.open(f)[0].data
            if np.percentile(d, 90)>1500 and np.percentile(d, 90)<40000:
                lfiles.append(f)

        if len(lfiles) == 0:
            print "WARNING!!! Could not find suitable flats for band %s"%b
            continue
        ffile ="lflat_"+b
        np.savetxt(ffile, np.array(lfiles), fmt="%s")
    
        
        #Cleaning of old files
        if(os.path.isfile(out)): os.remove(out)
        if(os.path.isfile(out_norm)): os.remove(out_norm)
        if(os.path.isfile("Flat_stats")): os.remove("Flat_stats")
        
        
        #Combine flats
        iraf.imcombine(input = "@"+ffile, \
                        output = out, \
                        combine = "median",\
                        scale = "mode",
                        weight = "exposure")
        iraf.imstat(out, fields="image,npix,mean,stddev,min,max,mode", Stdout="Flat_stats")
        st = np.genfromtxt("Flat_stats", names=True, dtype=None)
        #Normalize flats
        iraf.imarith(out, "/", st["MODE"], out_norm)
        
        #Do some cleaning
        print 'Removing from lfiles'
        for f in glob.glob('b_*_%s.fits'%b):
            os.remove(f)

        os.remove(ffile)
        
        
        if os.path.isfile(fsfile):
            os.remove(fsfile)
        if os.path.isfile(fffile):
            os.remove(fffile)

def reduce_image(img, flatdir=None, biasdir=None, cosmic=True, astrometry=True, channel='rc', target_dir='reduced'):
    '''
    Applies Flat field and bias calibrations to the image.
    
    Steps:
    
    1. - Solve astrometry on the entire image.
    2. - Compute master bias and de-bias the image.
    3. - Separate the image into 4 filters.
    4. - Compute flat field for each filter and apply flat fielding on the image.
    5. - Computes cosmic ray rejectionon the entire image.
    6. - Compute zeropoint for each image and store in a log file.
    7. - Plot zeropoint for the night.
    '''
    
    print "Reducing image ", img    

    objectname = fitsutils.get_par(img, "OBJECT").replace(" ", "").replace("]","").replace("[", "")

    print "For object", objectname
    
    #Change to image directory
    mydir = os.path.dirname(img)
    if mydir=="": mydir = "."
    mydir = os.path.abspath(mydir)
    os.chdir(mydir)
    #Create destination directory
    if (not os.path.isdir(target_dir)):
        os.makedirs(target_dir)

    #Rename to the image name only
    img = os.path.basename(img)


    if (astrometry):
        print "Solving astometry for the whole image..."
        img = solve_astrometry(img)
        astro = "a_"
    else:
        astro = ""
        
    
    #Compute BIAS
    if (biasdir == None or biasdir==""): biasdir = "."
    create_masterbias(biasdir)
    
    bias_slow = os.path.join(biasdir, "Bias_%s_%s.fits"%(channel, 'slow'))
    bias_fast = os.path.join(biasdir, "Bias_%s_%s.fits"%(channel, 'fast'))
    
    # Running IRAF to DE-BIAS
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    
    #Compute flat field
    if (flatdir == None or flatdir==""): flatdir = "."
    create_masterflat(flatdir, biasdir)
    
    #New names for the object.
    debiased = "b_" + astro + img
    print "Creating debiased file, ",debiased
    
    if (not os.path.isfile(bias_slow) or not os.path.isfile(bias_fast)):
        print "Master bias not found!"
        return

    #Debias
    if (fitsutils.get_par(img, "ADCSPEED")==2):
        iraf.imarith(img, "-", bias_fast, debiased)
        fitsutils.update_par(debiased, "BIASFILE", bias_fast)
        fitsutils.update_par(debiased, "RDNOISE", 20.)

    else:
        iraf.imarith(img, "-", bias_slow, debiased)
        fitsutils.update_par(debiased, "BIASFILE", bias_slow)
        fitsutils.update_par(debiased, "RDNOISE", 4.)

    #Set negative counts to zero
    hdu = pf.open(debiased)
    header = hdu[0].header
    hdu[0].data[hdu[0].data<0] = 0
    hdu.writeto(debiased, clobber=True)

    #Slicing the image for flats  
    slice_names = slice_rc(debiased)

    
    #Remove un-sliced image
    os.remove(debiased)

    # DE-flat each filter and store under object name
    for i, debiased_f in enumerate(slice_names):
        b = fitsutils.get_par(debiased_f, 'filter')
        
        deflatted = "f_b_" + astro + objectname + "_%s.fits"%b

        #Flat to be used for that filter
        flat = os.path.join(flatdir, "Flat_%s_%s_norm.fits"%(channel, b))

        if (not os.path.isfile(flat)):
            print "Master flat not found in", flat
            return
        #Cleans the deflatted file if exists
        if (os.path.isfile(deflatted)):
            os.remove(deflatted)
            
        iraf.imarith(debiased_f, "/", flat, deflatted)
        
        #Removes the de-biased file
        os.remove(debiased_f)
        
        print "Updating header with original filename and flat field used."
        fitsutils.update_par(deflatted, "ORIGFILE", img)
        fitsutils.update_par(deflatted, "FLATFILE", flat)

        slice_names[i] = deflatted
            
            
    if (cosmic):
        print "Correcting for cosmic rays..."
        # Correct for cosmics each filter
        for i, deflatted in enumerate(slice_names):
            cclean = "c_" +name
            clean_cosmic(os.path.join(os.path.abspath(mydir), deflatted), cclean)
            slice_names[i] = cclean
           
    #Moving files to the target directory
    for name in slice_names:
        if (os.path.isfile(name)):
            shutil.move(name, os.path.join(target_dir, name))
    
