# -*- coding: utf-8 -*-
"""
Created on Sun Jun 14 13:11:52 2015

@author: nadiablago
"""
import fitsutils
import os, glob, shutil
import numpy as np
try:
    from pyraf import iraf 
except:
    pass
import pyfits as pf
from matplotlib import pylab as plt
import subprocess
import argparse
import shutil
import time_utils
import time
import datetime
import logging
import sextractor

#Log into a file
FORMAT = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
#root_dir = "/tmp/logs/"
root_dir = "/scr2/sedm/logs/"
now = datetime.datetime.utcnow()
timestamp=datetime.datetime.isoformat(now)
timestamp=timestamp.split("T")[0]
logging.basicConfig(format=FORMAT, filename=os.path.join(root_dir, "rcred_{0}.log".format(timestamp)), level=logging.INFO)
logger = logging.getLogger('rcred')
    
    
def create_masterbias(biasdir=None, channel='rc'):
    '''
    Combines slow and fast readout mode biases for the specified channel.
    '''
    
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    
    if (biasdir is None) or biasdir=="": biasdir = "."
        
    outs = "Bias_%s_slow.fits"%channel
    outf = "Bias_%s_fast.fits"%channel

    doslow = True
    dofast = True
    if (os.path.isfile(os.path.join(biasdir,outs))): 
        logger.warn( "%s master Bias exists!"%outs)
        doslow = False
    if ( os.path.isfile(os.path.join(biasdir,outf))):
        logger.warn("%s master Bias exists!"%outs)
        dofast = False
    if(doslow or dofast):
        logger.info("Starting the Master Bias creation!")
    else:
        return

    os.chdir(biasdir)        
        
    lfastbias = []
    lslowbias = []
    
    #Select all filts that are Bias with same instrument
    for f in glob.glob("rc*fits"):
        try:
            if ( "BIAS" in str.upper(fitsutils.get_par(f, "IMGTYPE")) ):
                if (fitsutils.get_par(f, "ADCSPEED")==2):
                    lfastbias.append(f)
                else:
                    lslowbias.append(f)
        except:
            pass
                
    logger.info("Files for bias SLOW mode: %s"% lslowbias)
    logger.info( "Files for bias FAST mode: %s"% lfastbias)
    
    if len(lfastbias) > 0 and dofast:
        bfile_fast ="lbias_fast_"+channel
        np.savetxt(bfile_fast, np.array(lfastbias), fmt="%s")
        if (os.path.isfile("Bias_stats_fast")): os.remove("Bias_stats_fast")
        iraf.imstat("@"+bfile_fast, Stdout="Bias_stats_fast")
        
        st = np.genfromtxt("Bias_stats_fast", names=True, dtype=None)
        logger.info("%s"%st)
        
        iraf.imcombine(input = "@"+bfile_fast, \
                    output = outf, \
                    combine = "median",\
                    scale = "mode")
        os.remove(bfile_fast)
        
        #copy into the reference folder with current date
        newdir = os.path.join("../../refphot/", os.path.basename(os.path.abspath(biasdir)))
        if (not os.path.isdir(newdir)):
            os.makedirs(newdir)
        shutil.copy(outf, os.path.join(newdir, os.path.basename(outf)) )
    else:
        copy_ref_calib(biasdir, outf)


    if len(lslowbias) > 0 and doslow:

        bfile_slow ="lbias_slow_"+channel
        np.savetxt(bfile_slow, np.array(lslowbias), fmt="%s")
        if (os.path.isfile("Bias_stats_slow")): os.remove("Bias_stats_slow")
        iraf.imstat("@"+bfile_slow, Stdout="Bias_stats_slow")
        
        st = np.genfromtxt("Bias_stats_slow", names=True, dtype=None)
        logger.info("%s"%st)
        
        iraf.imcombine(input = "@"+bfile_slow, \
                    output = outs, \
                    combine = "median",\
                    scale = "mode")
        os.remove(bfile_slow)
        
        #copy into the reference folder with current date
        newdir = os.path.join("../../refphot/", os.path.basename(os.path.abspath(biasdir)))
        if (not os.path.isdir(newdir)):
            os.makedirs(newdir)
        shutil.copy(outs, os.path.join(newdir, os.path.basename(outs)) )  
    else:
        copy_ref_calib(biasdir, outs)

def create_masterflat(flatdir=None, biasdir=None, channel='rc'):
    '''
    Creates a masterflat from both dome flats and sky flats if the number of counts in the given filter
    is not saturated and not too low (between 3000 and 40000). 
    '''
    
    
    if (flatdir == None or flatdir==""): flatdir = "."
        
    if (biasdir == None or biasdir==""): biasdir = "."
        
    os.chdir(flatdir)
    
    if (len(glob.glob("Flat_%s*norm.fits"%channel)) == 4):
        logger.info( "Master Flat exists!")
        return 
    if (len(glob.glob("Flat_%s*norm.fits"%channel)) > 0):
        logger.info( "Some Master Flat exist!")
        return 
    else:
        logger.info( "Starting the Master Flat creation!")

    bias_slow = "Bias_%s_slow.fits"%channel
    bias_fast = "Bias_%s_fast.fits"%channel
    
    if (not os.path.isfile(bias_slow) and not os.path.isfile(bias_fast) ):
        create_masterbias(biasdir)
     
    lsflat = []
    lfflat = []
    
    obj = ""
    imtype = ""
    
    #Select all filts that are Flats with same instrument
    for f in glob.glob(channel+"*fits"):
        try:
            if fitsutils.has_par(f, "OBJECT"):
                obj = str.upper(fitsutils.get_par(f, "OBJECT"))
            else:
                continue

            if fitsutils.has_par(f, "IMGTYPE"):
                imtype = str.upper(fitsutils.get_par(f, "IMGTYPE"))
            else:
                continue
        
            #if ("RAINBOW CAM" in str.upper(fitsutils.get_par(f, "CAM_NAME")) and  ("DOME" in  obj or "FLAT" in obj or "Twilight" in obj or "TWILIGHT" in imtype or "DOME" in imtype)):
            if ( "TWILIGHT" in imtype):

                if (fitsutils.get_par(f, "ADCSPEED")==2):
                    lfflat.append(f)
                else:
                    lsflat.append(f)
        except:
            logger.error( "Error with retrieving parameters for file %s"% f)
            pass
                
    logger.info( "Files for slow flat %s"% lsflat)
    logger.info( "Files for fast flat %s"% lfflat)
    
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
    
    #Remove the list files
    os.remove(fsfile)
    os.remove(fffile)
    
    #Slices the flats.
    debiased_flats = glob.glob("b_*.fits")
    for f in debiased_flats:
        logger.info( "Slicing file %s"% f)
        slice_rc(f)
        #Remove the un-sliced file
        os.remove(f)
        
    #Selects the ones that are suitable given the number of counts and combines them.
    bands = ['u', 'g', 'r', 'i']
    for b in bands:
        out = "Flat_%s_%s.fits"%(channel, b)
        out_norm = out.replace(".fits","_norm.fits")

        if (os.path.isfile(out_norm)):
            logger.error( "Master Flat for filter %s exists. Skipping..."%b)
            continue
        
        lfiles = []
        for f in glob.glob('b_*_%s.fits'%b):
            d = pf.open(f)[0].data
            if np.percentile(d, 90)>4000 and np.percentile(d, 90)<40000:
                lfiles.append(f)

        if len(lfiles) == 0:
            logger.error( "WARNING!!! Could not find suitable flats for band %s"%b)
            continue
        if len(lfiles) < 3:
            logger.error( "WARNING!!! Could find less than 3 flats for band %s. Skipping, as it is not reliable..."%b)
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
        logger.info( 'Removing from lfiles')
        for f in glob.glob('b_*_%s.fits'%b):
            os.remove(f)

        os.remove(ffile)
        
        
        if os.path.isfile(fsfile):
            os.remove(fsfile)
        if os.path.isfile(fffile):
            os.remove(fffile)

        #copy into the reference folder with current date
        newdir = os.path.join("../../refphot/", os.path.basename(os.path.abspath(flatdir)))
        if (not os.path.isdir(newdir)):
            os.makedirs(newdir)
        shutil.copy(out_norm, os.path.join(newdir, os.path.basename(out_norm)) )   


def get_median_bkg(img):
    '''
    Computes the median background.
    '''
    hdu = pf.open(img)
    header = hdu[0].header
    bkg = np.median(hdu[0].data[hdu[0].data > 0])
    return bkg
    
def solve_edges(flat, band):
    '''
    Finds the edges of the flat, and extends the respons of the average to the edges.
    '''
    
    edgedic = {"u":[], "g":[], "i":[], "r":[]}
    
    
    
def copy_ref_calib(curdir, calib="Flat"):
    '''
    Reference master Bias and master Flat are stored in the refphot folder.
    The files are copied if they are not found in the folder where the photometry is being reduced.
    '''

    #Check which dates have the calibration files that we need.
    listcalib = glob.glob(os.path.join(curdir, "../../refphot/*/", calib+"*"))
    #Obtain the folder name
    listdirs = [os.path.dirname(l) for l in listcalib]    
    #Unique name
    calibdates = np.array(list(set(listdirs)))
    
    #Get the date of the current directory
    curdir = os.path.abspath(curdir)
    
    #compile the dates in datetime format
    dates = []
    for f in calibdates:
        dateobs = os.path.basename(f)
        dates.append(datetime.datetime.strptime(dateobs, "%Y%m%d"))

    dates = np.array(dates)    
    curdate = datetime.datetime.strptime(os.path.basename(curdir), "%Y%m%d")
    
    #Select the folder that is closes to the date of the current directory
    lastdir = np.array(calibdates)[np.argmin(np.abs(curdate-dates))]
    
    #Copy all calibration files that match the calib filter.
    for c in glob.glob(os.path.join(lastdir, calib+"*")):    
        shutil.copy(c, os.path.join(curdir, os.path.basename(c)))     
        
        

def solve_astrometry(img, radius=3, with_pix=True, overwrite=False, tweak=3):
    '''
    img: fits image where astrometry should be solved.
    radius: radius of uncertainty on astrometric position in image.
    '''

    ra = fitsutils.get_par(img, 'OBJRA')
    dec = fitsutils.get_par(img, 'OBJDEC')
    logger.info( "Solving astrometry on field with (ra,dec)=%s %s"%(ra, dec))
    
    astro = "a_" + img
    
    #If astrometry exists, we don't run it again.
    if (os.path.isfile(astro) and not overwrite):
        return astro
        
    #Store a temporary file with the multiplication of the mask
    mask = "mask.fits"

    cmd = "solve-field --ra %s --dec %s --radius %.4f -p --new-fits %s \
      -W none -B none -P none -M none -R none -S none -t %d --overwrite %s "%(ra, dec, radius, astro, tweak, img)
    if (with_pix):
        cmd = cmd + " --scale-units arcsecperpix  --scale-low 0.375 --scale-high 0.4"
    logger.info( cmd)

    subprocess.call(cmd, shell=True)
    
    #Cleaning after astrometry.net
    if (os.path.isfile(img.replace(".fits", ".axy"))):
        os.remove(img.replace(".fits", ".axy"))
    if (os.path.isfile(img.replace(".fits", "-indx.xyls"))):
        os.remove(img.replace(".fits", "-indx.xyls"))
    if (os.path.isfile("none")):
        os.remove("none")
        
    is_on_target(img)
    
    return astro
    

def make_mask_cross(img):  
    
    maskname = os.path.join(os.path.dirname(img), "mask.fits")
    
    if (os.path.isfile(maskname)):
        return maskname
        
    f = pf.open(img)
    data = f[0].data
    
    corners = {
    "g" : [1, 850, 1, 850],
    "i" : [1, 850, 110, 2045],
    "r" : [1150, 2045, 1150, 2045],
    "u" : [1130, 2045, 1, 850]
    }
    
    newdata = np.ones_like(data)*np.nan
    
    for band in corners.keys():
        ix = corners[band]
        newdata[ix[0]:ix[1], ix[2]: ix[3]] = 1
        
    f[0].data = newdata
    f.writeto(maskname, clobber=True)
    
    return maskname
    
def get_masked_image(img):

    # Running IRAF
    iraf.noao(_doprint=0)  
    mask = make_mask_cross(img)    
    masked = img.replace(".fits", "_masked.fits")
    
    if (os.path.isfile(masked)):
        os.remove(masked)
    iraf.imarith(img, "*", mask, masked)
    
    return masked

    
def slice_rc(img):
    '''
    Slices the Rainbow Camera into 4 different images and adds the 'filter' keyword in the fits file.
    '''
    fname = os.path.basename(img)
    fdir = os.path.dirname(img)    
    
    # Running IRAF
    iraf.noao(_doprint=0)    
    
    corners = {
    "g" : [1, 910, 1, 900],
    "i" : [1, 910, 1060, 2045],
    "r" : [1040, 2045, 1015, 2045],
    "u" : [1030, 2045, 1, 900]
    }
    
    
    filenames = []
        
    for i, b in enumerate(corners.keys()):
        logger.info( "Slicing for filter %s"% b)
        name = fname.replace(".fits", "_%s.fits"%b)
        
        #Clean first
        if (os.path.isfile(name)):
            os.remove(name)
            
        iraf.imcopy("%s[%d:%d,%d:%d]"%(img, corners[b][0], corners[b][1], corners[b][2], corners[b][3]), name)
         
        fitsutils.update_par(name, 'filter', b)
        is_on_target(name)
        
        filenames.append(name)
    
    return filenames

def is_on_target(image):
    '''
    Add as a parameter whether the image is on target or not.
    
    '''
    import pywcs
    import coordinates_conversor as cc
    
    ra, dec = cc.hour2deg(fitsutils.get_par(image, 'OBJRA'), fitsutils.get_par(image, 'OBJDEC'))

    impf = pf.open(image)
    wcs = pywcs.WCS(impf[0].header)
    #pra, pdec = wcs.wcs_sky2pix(ra, dec, 1)
    pra, pdec = wcs.wcs_sky2pix(np.array([ra, dec], ndmin=2), 1)[0]

    shape = impf[0].data.shape
    
    if (pra > 0)  and (pra < shape[0]) and (pdec > 0) and (pdec < shape[1]):
        fitsutils.update_par(image, "ONTARGET", 1)
    else:
        fitsutils.update_par(image, "ONTARGET", 0)
        
    
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
    

    
                
def get_overscan_bias_rc(img):
    '''
    Bias from overscan region.
    '''
    f = pf.open(img)
    bias = np.nanmedian(f[0].data[990-100:990+100,970-100:970+100].flatten())
    
    return bias
    
def get_sequential_name(target_dir, name, i=0):
        '''
        Gets a sequential name if we have a file with the same object name imaged several times.
        '''
        newname = os.path.join(target_dir, name).replace(".fits", "_%d.fits"%i)
        #If the destination file does not exist...
        if (os.path.isfile(name) and not os.path.isfile(newname)):
            return newname
        #If it does exist, but it is another exposure
        elif(os.path.isfile(name) and os.path.isfile(newname)):
            #If it is not the same exposure, add with different name. Otherwise, replace.
            if fitsutils.get_par(name, "JD") != fitsutils.get_par(newname, "JD"):
                newname = get_sequential_name(target_dir, name, i=i+1)
        
        return newname      

def init_header_reduced(image):
    '''
    IQWCS = 1 or 0 / Indicates astrometry has been solved for the field
    IQZEROPT =  1 or 0 /indicates if the zero point was calculated for the image    
    SKYBKG = FLOAT / Average sky background given in counts    
    SEEPIX = FLOAT  / Seeing expressed in pixels    
    ZPCAT = 'String' / Catalog used to calculate zero point    
    ZEROPTU = FLOAT  / Zero point uncertainty    
    ZEROPT  = 'FLOAT'  / Zero point by comparison to catalog
    '''
    
    pardic = {"IQWCS" : 0,\
                "IQZEROPT" : 0,\
                "SKYBKG": 0,\
                "SEEPIX": 0,\
                "ZPCAT" : "none",\
                "ZEROPTU" : 0.,\
                "ZEROPT" : 0.}
    fitsutils.update_pars(image, pardic)
    
    
def reduce_image(image, flatdir=None, biasdir=None, cosmic=False, astrometry=True, channel='rc', target_dir='reduced', overwrite=False):
    '''
    Applies Flat field and bias calibrations to the image.
    
    Steps:
    
    1. - Solve astrometry on the entire image.
    2. - Compute master bias (if it does not exist) and de-bias the image.
    3. - Separate the image into 4 filters.
    4. - Compute flat field for each filter (if it does not exist) and apply flat fielding on the image.
    5. - Computes cosmic ray rejectionon the entire image.

    '''
    
    logger.info("Reducing image %s"% image)    

    print "Reducing image ", image    

    #Initialize the basic parameters.
    init_header_reduced(image)
    
    
    
    imname = os.path.basename(image).replace(".fits", "")
    try:
        objectname = fitsutils.get_par(image, "NAME").replace(" ","")+"_"+fitsutils.get_par(image, "FILTER")
    except:
        logger.error( "ERROR, image "+ image + " does not have a NAME or a FILTER!!!")
        return

    print "For object", objectname
    logger.info( "For object %s"% objectname)

        
    #Get basic statistics for the image
    nsrc, fwhm, ellip = sextractor.get_image_pars(image)
    
    logger.info( "Sextractor statistics: nscr %d, fwhm (pixel) %.2f, ellipticity %.2f"% (nsrc, fwhm, ellip))
    print "Sextractor statistics: nscr %d, fwhm (pixel) %.2f, ellipticity %.2f"% (nsrc, fwhm, ellip)

    
    dic = {"SEEPIX": fwhm/0.394, "NSRC":nsrc, "ELLIPTICITY":ellip}
    #Update the seeing information from sextractor
    fitsutils.update_pars(image, dic)

    
    
    #Change to image directory
    mydir = os.path.dirname(image)
    if mydir=="": mydir = "."
    mydir = os.path.abspath(mydir)
    os.chdir(mydir)
    #Create destination directory
    if (not os.path.isdir(target_dir)):
        os.makedirs(target_dir)

    #If we don't want to overwrite the already extracted images, we check wether they exist.
    if (not overwrite):
        existing = True
        for band in ['u', 'g', 'r', 'i']:
            destfile = os.path.join(target_dir, imname + "_f_b_a_%s_%s_0.fits"%(objectname, band))
            logger.info( "Looking if file %s exists"%( destfile, \
                (os.path.isfile(destfile) ) ) )
            existing = existing and (os.path.isfile(os.path.join(target_dir, "f_b_a_%s_%s_0.fits"%(objectname, band))) )
        if existing:
            return []


    #Rename to the image name only
    image = os.path.basename(image)


    astro = ""
    if (astrometry):
        logger.info( "Solving astometry for the whole image...")
        img = solve_astrometry(image)
        if (os.path.isfile(img)):
            astro="a_"
            fitsutils.update_par(img, "IQWCS", 1)
        else:
            logger.error( "ASTROMETRY DID NOT SOLVE ON IMAGE %s"% image)
            img = image


        
    
    #Compute BIAS
    if (biasdir is None or biasdir==""): biasdir = "."
    create_masterbias(biasdir)
    
    bias_slow = os.path.join(biasdir, "Bias_%s_%s.fits"%(channel, 'slow'))
    bias_fast = os.path.join(biasdir, "Bias_%s_%s.fits"%(channel, 'fast'))
    
    # Running IRAF to DE-BIAS
    iraf.noao(_doprint=0)
    iraf.imred(_doprint=0)
    iraf.ccdred(_doprint=0)
    
    #Compute flat field
    if (flatdir is None or flatdir==""): flatdir = "."
    create_masterflat(flatdir, biasdir)
    
    #New names for the object.
    debiased = "b_" + img
    logger.info( "Creating debiased file, %s"%debiased)
    
    if ( (fitsutils.get_par(img, "ADCSPEED")==0.1 and not os.path.isfile(bias_slow)) \
        or (fitsutils.get_par(img, "ADCSPEED")==2 and not os.path.isfile(bias_fast)) ):
        logger.warn( "Master bias not found! Tryting to copy from reference folder...")
        copy_ref_calib(mydir, "Bias")
        if ( (fitsutils.get_par(img, "ADCSPEED")==0.1 and not os.path.isfile(bias_slow)) \
        or (fitsutils.get_par(img, "ADCSPEED")==2 and not os.path.isfile(bias_fast)) ):
            logger.error( "Bias not found in reference folder")
            return


    #Clean first
    if (os.path.isfile(debiased)):
        os.remove(debiased)
        
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
    print "Creating sliced files, ", slice_names

    
    #Remove un-sliced image
    os.remove(debiased)

    # DE-flat each filter and store under object name
    for i, debiased_f in enumerate(slice_names):
        b = fitsutils.get_par(debiased_f, 'filter')
        
        deflatted = imname + "_f_b_" + astro + objectname + "_%s.fits"%b

        #Flat to be used for that filter
        flat = os.path.join(flatdir, "Flat_%s_%s_norm.fits"%(channel, b))

        if (not os.path.isfile(flat)):
            logger.warn( "Master flat not found in %s"% flat)
            copy_ref_calib(mydir, "Flat_%s_%s_norm"%(channel, b))
            continue
        else:
            logger.info( "Using flat %s"%flat)
            
        #Cleans the deflatted file if exists
        if (os.path.isfile(deflatted)):
            os.remove(deflatted)

        if (os.path.isfile(debiased_f) and os.path.isfile(flat)):
            logger.info( "Storing de-flatted %s as %s"%(debiased_f, deflatted))
            time.sleep(1)
            iraf.imarith(debiased_f, "/", flat, deflatted)
        else:
            logger.error( "SOMETHING IS WRONG. Error when dividing %s by the flat field %s!"%(debiased_f, flat))
        
        #Removes the de-biased file
        os.remove(debiased_f)
        
        logger.info( "Updating header with original filename and flat field used.")
        fitsutils.update_par(deflatted, "ORIGFILE", img)
        fitsutils.update_par(deflatted, "FLATFILE", flat)

        slice_names[i] = deflatted
            
            
    if (cosmic):
        logger.info( "Correcting for cosmic rays...")
        # Correct for cosmics each filter
        for i, deflatted in enumerate(slice_names):
            cclean = "c_" +img
            clean_cosmic(os.path.join(os.path.abspath(mydir), deflatted), cclean)
            slice_names[i] = cclean
                
            
    reduced_imgs = []   
    #Moving files to the target directory
    for name in slice_names:
        newname = get_sequential_name(target_dir, name, 0)
        bkg = get_median_bkg(name)
        fitsutils.update_par(name, "SKYBKG", bkg)
        shutil.move(name, newname)
        reduced_imgs.append(newname)
        
    return reduced_imgs

                
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

        Reduces the photometric images from SEDM Rainbow Camera.
        Requires either a list of images to be reduced, or a directory name where the night photometry is.
        As an option, can run lacosmic to remove cosmic rays.
        By default it invokes astrometry.net before the reduction.
        
        %run rcred.py -d PHOTDIR 
                
        Reduced images are stored in directory called "reduced", within the main directory.
        
        Optionally, it can be used to clean the reduction products generated by this pipeline within PHOTDIR diretory using -c option (clean).
        
        %run rcred.py -d PHOTDIR -c
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-l', '--filelist', type=str, help='File containing the list of fits for the night.',default=None)
    parser.add_argument('-d', '--photdir', type=str, help='Directory containing the science fits for the night.', default=None)
    parser.add_argument('-c', '--clean', action="store_true", help='Clean the reduced images?', default=False)
    parser.add_argument('-o', '--overwrite', action="store_true", help='re-reduce and overwrite the reduced images?', default=False)
    parser.add_argument('--cosmic', action="store_true", default=False, help='Whether cosmic rays should be removed.')

    args = parser.parse_args()
    
    filelist = args.filelist
    photdir = args.photdir
    cosmic = args.cosmic
    clean =  args.clean
    overwrite = args.overwrite
    
    myfiles = []


    if (not photdir is None and clean):
        for f in glob.glob("Flat*"):
            os.remove(f)
        for f in glob.glob("Bias*"):
            os.remove(f)
        for f in glob.glob("a_*fits"):
            os.remove(f)
        if (os.path.isdir(os.path.join(photdir, "reduced"))):
            shutil.rmtree(os.path.join(photdir, "reduced"))

    elif (not filelist is None ):
        mydir = os.path.dirname(filelist)
        if (mydir==""):
            mydir = "."
        os.chdir(mydir)
        
        myfiles = np.genfromtxt(filelist, dtype=None)
        
    elif (not photdir is None ):
        mydir = os.path.abspath(photdir)
        #Gather all RC fits files in the folder with the keyword IMGTYPE=SCIENCE
        for f in glob.glob(os.path.join(mydir, "rc*fits")):
            if (fitsutils.has_par(f, "IMGTYPE") and fitsutils.get_par(f, "IMGTYPE")=="SCIENCE"):
                myfiles.append(f)
    else:
        print '''ERROR! You should specify one of the following:\n
                   - A filelist name with the images you want to reduce [-l]
                   OR
                   - The name of the directory which you want to reduce [-d].'''

    if (False): #len(myfiles)>0):
    	create_masterbias()
    	create_masterflat()  
    
    #Reduce them
    for f in myfiles:
        print f
        make_mask_cross(f)
        if(fitsutils.has_par(f, "IMGTYPE") and fitsutils.get_par(f, "IMGTYPE") == "SCIENCE" or fitsutils.get_par(f, "IMGTYPE") == "ACQUISITION"):
	    try:
            	reduced = reduce_image(f, cosmic=cosmic, overwrite=overwrite) 
	    except:
		pass
    
