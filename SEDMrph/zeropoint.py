# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 14:57:05 2015

@author: nadiablago
"""
import matplotlib
matplotlib.use("Agg")
import zscale
import sys, urllib
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
import math, glob
import pywcs 
from astropy.io.votable import parse_single_table
import coordinates_conversor
import app_phot
from numpy.lib import recfunctions as rfn
import scipy.optimize as opt
import fitsutils
import os, shutil
import transformations
from astropy import stats
import argparse
import logging
import datetime


from ConfigParser import SafeConfigParser
import codecs

parser = SafeConfigParser()

configfile = os.environ["SEDMCONFIG"]

# Open the file with the correct encoding
with codecs.open(configfile, 'r') as f:
    parser.readfp(f)

_logpath = parser.get('paths', 'logpath')
_photpath = parser.get('paths', 'photpath')


FORMAT = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
now = datetime.datetime.utcnow()
timestamp=datetime.datetime.isoformat(now)
creationdate = timestamp
timestamp=timestamp.split("T")[0]

try:
    #Log into a file
    root_dir = _logpath
    logging.basicConfig(format=FORMAT, filename=os.path.join(root_dir, "rcred_{0}.log".format(timestamp)), level=logging.INFO)
    logger = logging.getLogger('zeropoint')
except:
    logging.basicConfig(format=FORMAT, filename=os.path.join("/tmp", "rcred_{0}.log".format(timestamp)), level=logging.INFO)
    logger= logging.getLogger("zeropoint")
    




def are_isolated(rav, decv, r):
    '''
    Returns a mask saying whether the stars with coordinates rav, decv are isolated in a radius of r arcsec.
    '''
    
    index = np.arange(len(rav))
    mask = []
    
    for i in np.arange(len(rav)):
        d = coordinates_conversor.get_distance(rav[i], decv[i], rav[index!=i], decv[index!=i])
        mask.append(~ np.any(d*3600<r))
    
    return np.array(mask)

def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    '''
    Produces a 2D gaussian centered in xo, yo with the parameters specified.
    xdata_tuple: coordinates of the points where the 2D Gaussian is computed.
    
    '''
    (x, y) = xdata_tuple                                                        
    xo = float(xo)                                                              
    yo = float(yo)                                                              
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)   
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)    
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)   
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))                                   
    return g.ravel()
    
def twoD_Gauss_test(theta=0):
    '''
    Generates a test Gaussian and fits it using the scisoft optimization software.
    '''
    # Create x and y indices
    x = np.linspace(0, 200, 201)
    y = np.linspace(0, 200, 201)
    x, y = np.meshgrid(x, y)
    
    #create data
    data = twoD_Gaussian((x, y), 3, 100, 100, 20, 40, theta, 10)
    
    # plot twoD_Gaussian data generated above
    plt.figure()
    plt.imshow(data.reshape(201, 201), origin="bottom", extent=(x.min(), x.max(), y.min(), y.max()))
    plt.colorbar()
    
    # add some noise to the data and try to fit the data generated beforehand
    initial_guess = (3,100,100,20,40,0,10)
    
    data_noisy = data + 0.2*np.random.normal(size=data.shape)
    
    popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), data_noisy, p0=initial_guess)
    
    data_fitted = twoD_Gaussian((x, y), *popt)
    
    fig, ax = plt.subplots(1, 1)
    ax.hold(True)
    ax.imshow(data_noisy.reshape(201, 201), cmap=plt.cm.jet, origin='bottom',
        extent=(x.min(), x.max(), y.min(), y.max()))
    ax.contour(x, y, data_fitted.reshape(201, 201), 8, colors='w')
    plt.show()

def find_fwhm(imfile, xpos, ypos, plot=True):
    '''
    Finds and returns the best parameters for the FWHM in arcsec for the stars marked with X, Y
    '''
    
    f = pf.open(imfile)
    img = f[0].data
    pix2ang = 0.394
    def_fwhm = 2./pix2ang
    rad = math.ceil(30./pix2ang)
    
    out = np.zeros(len(xpos), dtype=[('detected', np.bool), ('fwhm', np.float), ('e', np.float)])
    
    for i, (x_i,y_i) in enumerate(zip(xpos, ypos)):

        x_i = int(x_i)
	y_i = int(y_i)
        hrad = int(math.ceil(rad/2.))
        sub = img[x_i-hrad:x_i+hrad, y_i-hrad:y_i+hrad]
        x = np.linspace(0, len(sub), len(sub))
        y = np.linspace(0, len(sub), len(sub))
        X, Y = np.meshgrid(x, y)
    
        #(xdata_tuple, amplitude, xo, yo, def_fwhm, def_fwhm, theta, offset):
        def_x = np.argmax(np.sum(sub, axis=0))
        def_y = np.argmax(np.sum(sub, axis=1))

        initial_guess = (100, def_x, def_y, def_fwhm, def_fwhm, 0, np.percentile(sub, 50))
        detected = True
        try:
            popt, pcov = opt.curve_fit(twoD_Gaussian, (X, Y), sub.flatten(), p0=initial_guess)
            fwhm_x = np.abs(popt[3])*2*np.sqrt(2*np.log(2))
            fwhm_y = np.abs(popt[4])*2*np.sqrt(2*np.log(2))
            amplitude=popt[0]
            background=popt[-1]
        except:
            detected = False
            fwhm_x = 0
            fwhm_y = 0
            amplitude = 0
            background = 0
        
        
        detected = detected * (amplitude> 0)
        
        if (not detected):
            fwhm_x = 0
            fwhm_y = 0     
        
        
        logger.info("%s %s Amplitude %.3f\t BG %.3f\t BG_stats %.3f\t  FWHM_x,FWHM_y=(%.3f, %.3f)"%(i, detected, amplitude, background, np.percentile(sub, 50), fwhm_x, fwhm_y))
        
        out[i] = (detected, np.average([fwhm_x, fwhm_y]), np.minimum(fwhm_x, fwhm_y) / np.maximum(fwhm_x, fwhm_y))
        
        if (detected & plot):
            data_fitted = twoD_Gaussian((X, Y), *popt)
            
            fig, (ax, ax2) = plt.subplots(1, 2)
            ax.hold(True)
            ax.imshow(sub, cmap=plt.cm.jet, origin='bottom', extent=(x.min(), x.max(), y.min(), y.max()))
            ax.contour(X, Y, data_fitted.reshape(sub.shape[0], sub.shape[1]), 5, colors='w')
            ax2.imshow(sub-data_fitted.reshape(sub.shape[0], sub.shape[1]), cmap=plt.cm.jet, origin='bottom', extent=(x.min(), x.max(), y.min(), y.max()))
            ax2.contour(X, Y, data_fitted.reshape(sub.shape[0], sub.shape[1]), 5, colors='w')
            ax.scatter(def_x, def_y, marker="*", s=100, color="yellow")
            plt.title("DETECTED for Position X,Y = %d,%d"%(x_i,y_i))
            plt.savefig(os.path.join(os.path.dirname(imfile), "gauss_%d"%i))
            plt.clf()
        if ((not detected) & plot):           
            fig, ax = plt.subplots(1)
            ax.hold(True)
            ax.imshow(sub, cmap=plt.cm.jet, origin='bottom', extent=(x.min(), x.max(), y.min(), y.max()))
            plt.title("NOT DETECTED for Position X,Y = %d,%d"%(x_i,y_i))
            plt.savefig(os.path.join(os.path.dirname(imfile), "gauss_%d"%i))
            plt.clf()
            #plt.show()
            
    return out
    
def clean_tmp_files():
    
    files = ["/tmp/tmp_sdss_%s.cat"%creationdate, '/tmp/tmp_%s.cat'%creationdate, '/tmp/tmp_apass_%s.cat'%creationdate, '/tmp/tmp_sdss_%s.cat'%creationdate, '/tmp/sdss_cat_det_%s.txt'%creationdate]
    
    for f in files:
        if os.path.isfile(f):
            os.remove(f)
    
def extract_star_sequence(imfile, band, plot=True, survey='sdss', debug=False, refstars=None, plotdir="."):
    '''
    Given a fits image: imfile and a the name of the band which we want to extract the sources from,
    it saves the extracted sources into  '/tmp/sdss_cat_det.txt' file.
    If the band does not match the bands in the survey, a change is performed to adapt to the new band.
    
    If plotting activated, plots the USNOB1 field of stars on top of the star field image.
    Red circles are stars identified from the catalogue in the correct magnitude range.
    Yellow circles are stars that are isolated.
    
    '''
    
    survey = str.lower(survey)
    minmag = 15
    maxmag = 21.0
        
    f = pf.open(imfile)
    wcs = pywcs.WCS(f[0].header)

    img = f[0].data
    img[img<0] = 0
    
    ra, dec = wcs.wcs_pix2sky(np.array([img.shape[0]/2, img.shape[1]/2], ndmin=2), 1)[0]
    ra0, dec0 = wcs.wcs_pix2sky(np.array([img.shape[0], img.shape[1]], ndmin=2), 1)[0]

    sr = 2*np.abs(dec-dec0)
    logger.info("%.4f %.4f %.4f"%( ra,dec, sr))
    
    if not refstars is None:
        shutil.copy(refstars, "/tmp/tmp_sdss_%s.cat"%creationdate)
        catalog = np.genfromtxt("/tmp/tmp_sdss_%s.cat"%creationdate, names=True, dtype=None, delimiter=",")
        cat_ra = catalog["ra"]
        cat_dec = catalog["dec"]
        try:
            mag = catalog["R"]
        except:
            mag = catalog["r"]
        
    elif str.lower(survey) =='usnob1':
        #ra, dec = coordinates_conversor.hour2deg(f[0].header['RA'], f[0].header['DEC'])
        #SEDM FoV is 6.5 arcmin, due to uncertainties in the position, 4 arcmin radius assumed.
        # Download USNO-B1 catalog for the position
        catalog_url = 'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=USNO-B1&RA=%.5f&DEC=%.5f&SR=%.4f&VERB=1' % (ra, dec, sr)
        logger.info( "Downloading USNO-B1 catalog...")
        urllib.urlretrieve(catalog_url, '/tmp/tmp_%s.cat'%creationdate)
        
        # Read RA, Dec and magnitude from XML format USNO catalog
        catalog = parse_single_table('/tmp/tmp_%s.cat'%creationdate)
        cat_ra = catalog.array['RA'].data
        cat_dec = catalog.array['DEC'].data
        cat_R1mag = catalog.array['R1'].data
        cat_R2mag = catalog.array['R2'].data
        cat_B1mag = catalog.array['B1'].data
        cat_B2mag = catalog.array['B2'].data
        cat_I2mag = catalog.array['I2'].data
        
        cat_R1mag[cat_R1mag==0] = np.nan
        cat_R2mag[cat_R2mag==0] = np.nan
        cat_B1mag[cat_B1mag==0] = np.nan
        cat_B2mag[cat_B2mag==0] = np.nan
        cat_I2mag[cat_I2mag==0] = np.nan
        
        Bmag = np.nanmean(np.array([cat_B1mag, cat_B2mag]), axis=0)    
        Rmag = np.nanmean(np.array([cat_R1mag, cat_R2mag]), axis=0)    
        Imag = cat_I2mag
            
        mag = Rmag
    elif survey == "apass":
        # Download USNO-B1 catalog for the position
        catalog_url = 'https://www.aavso.org/cgi-bin/apass_download.pl?ra=%.5f&dec=%.5f&radius=%.4f8&outtype=1' % (ra, dec, sr)
        print "Downloading APASS catalog..."
        urllib.urlretrieve(catalog_url, '/tmp/tmp_apass_%s.cat'%creationdate)
        catalog = np.genfromtxt('/tmp/tmp_apass_%s.cat'%creationdate, delimiter=",", names=True)
        if (np.ndim(catalog)==0):
		return False
        cat_ra = catalog['radeg']
        cat_dec = catalog['decdeg']
        mag = catalog['Sloan_r']
        
    elif (survey=='sdss'):
        minmag = 15
        maxmag = 21.5
        catalog_url='http://skyserver.sdss.org/dr9/en/tools/search/x_radial.asp?ra=%.5f&dec=%.5f&check_type=type&type=6&radius=%.4f&check_u=u&min_u=%.2f&max_u=%.2f&check_g=g&min_g=%.2f&max_g=%.2f&check_r=r&min_r=%.2f&max_r=%.2f&check_i=i&min_i=%.2f&max_i=%.2f&check_z=z&min_z=%.2f&max_z=%.2f&entries=top&topnum=500&format=csv'%\
            (ra, dec, sr*60,minmag,maxmag,minmag,maxmag,minmag,maxmag,minmag,maxmag,minmag,maxmag)
        logger.info( "Downloading SDSS catalog...")
        logger.info( "%s"%catalog_url )
        urllib.urlretrieve(catalog_url, '/tmp/tmp_sdss_%s.cat'%creationdate)
        catalog = np.genfromtxt("/tmp/tmp_sdss_%s.cat"%creationdate, delimiter=",", names=True)

	if (np.ndim(catalog)==0):
	    return False
        try:
            cat_ra = np.array(catalog['ra'], ndmin=1)
            cat_dec = np.array(catalog['dec'], ndmin=1)
            if (band in catalog.dtype.names):
                print "SDSS filter detected"
                mag = np.array(catalog[band], ndmin=1)
            elif(band in ["U", "B", "V", "R", "I", "Z"]):
                print "Johnson filter detected."
                john = transformations.sdss2johnson('/tmp/tmp_sdss_%s.cat'%creationdate, savefile='/tmp/tmp_sdss_%s.cat'%creationdate)
                mag = john[band]
                catalog = np.genfromtxt('/tmp/tmp_sdss_%s.cat'%creationdate, dtype=None, names=True, delimiter=",")
            else:
                print "Unknown band!!", band
        except IOError:
            logger.error( "Problems with SDSS image %s"% band)
            return False
        except ValueError:
            logger.error( "Problems with the catalogue for the image")
            return False

    #Convert ra, dec position of all stars to pixels.
    star_pix = np.array([0,0])
    for i in range(len(cat_ra)):
        # Get pixel coordinates of USNO stars
        s = wcs.wcs_sky2pix(np.array([cat_ra[i], cat_dec[i]], ndmin=2), 1)[0]
        star_pix = np.row_stack((star_pix, s))
    
    star_pix = star_pix[1:]
    pix2ang = 0.394
    rad = math.ceil(25./pix2ang)
    #Select only the stars within the image.
    mask = (star_pix[:,0]>-rad) * (star_pix[:,0]<img.shape[1]+rad)*(star_pix[:,1]>-rad) * (star_pix[:,1]<img.shape[0]+rad)
    if (band == 'u'):    
        mask = mask * (mag < 19)
    #Select only stars isolated in a radius of ~12 arcsec.
    mask2 = np.array(are_isolated(cat_ra[mask], cat_dec[mask], 15.))
    if (len(mask2)==0):
	logger.error("No good stars left")
	return False   
 
 
    #Select only stars that are within the proper magnitude range
    mask3 = (mag[mask][mask2] < maxmag) * (mag[mask][mask2] > minmag) 
    
    mask3 = mask3 * (star_pix[:,0][mask][mask2]>rad) * (star_pix[:,0][mask][mask2]<img.shape[1]-rad)*(star_pix[:,1][mask][mask2]>rad) * (star_pix[:,1][mask][mask2]<img.shape[0]-rad)


    if (survey=='usnob1'):
        output = np.column_stack((star_pix[:,0][mask][mask2][mask3], star_pix[:,0][mask][mask2][mask3], \
        cat_ra[mask][mask2][mask3], cat_dec[mask][mask2][mask3], Bmag[mask][mask2][mask3], Rmag[mask][mask2][mask3], Imag[mask][mask2][mask3]))
        np.savetxt('/tmp/usnob1_cat.txt', output, fmt="%.5f %.5f %.5f %.5f %.5f %.5f %.5f", header='#X Y ra dec B R I')
        print "Saved to", '/tmp/usnob1_cat.txt'
    elif (survey=='apass'):
            
        if not np.any(mask2) and not np.any(mask3):
            print star_pix
            print "No stars left...", mask, mask2, mask3
            return
        else:
            print catalog[mask][mask2][mask3]

        catalog = catalog[mask][mask2][mask3]
        s = star_pix[mask][mask2][mask3]
        z = np.zeros(len(s), dtype=[('x','f8'), ('y', 'f8')])
        z['x'] = s[:,0]
        z['y'] = s[:,1]
                
        for n  in catalog.dtype.names:
            z = rfn.append_fields(z, names=n, data=catalog[n], usemask=False)
           
        fmt = "%.5f"
        for i in range(len(z[0])-1):
            fmt +=  " %.5f"

        np.savetxt('/tmp/apass_cat.txt', z, fmt=fmt, \
        header='x y radeg raerr decdeg decerr number_of_Obs V dV B dB g dg r dr i di')
        print "Saved to", '/tmp/apass_cat.txt'
    elif survey=='sdss':
        if (not np.any(mask) and not np.any(mask2) and not np.any(mask3)) or len(catalog[mask][mask2][mask3])==0:
            print star_pix
            print "No stars left...", mask, mask2, mask3
            return False
        else:
            catalog = catalog[mask][mask2][mask3]
            s = star_pix[mask][mask2][mask3]

            print "left %d stars"%(len(catalog)), catalog.dtype.names

            z = np.zeros(len(s), dtype=[('x','f8'), ('y', 'f8')])
            z['x'] = s[:,0]
            z['y'] = s[:,1]
                    
            header='x y '
            for n  in catalog.dtype.names:
                if n in ["objid", "ra", "dec", "u", "g", "r", "i", "z", "Err_u", "Err_g", "Err_r", "Err_i", "Err_z"] or n in ['id', 'ra', 'dec', 'U', 'B', 'V', 'R', 'I', 'dU', 'dB', 'dV', 'dR', 'dI']:
                    z = rfn.append_fields(z, names=n, data=catalog[n], usemask=False)
                    header += n.replace("Err_", "d") + " "
               
            fmt = "%.5f"
            for i in range(len(z[0])-1):
                fmt +=  " %.5f"
    
            np.savetxt('/tmp/sdss_cat_%s.txt'%creationdate, z, fmt=fmt, header = header)
            logger.info( "Saved catalogue stars to %s"% ('/tmp/sdss_cat_%s.txt'%creationdate))
            
            #Find FWHM for this image            
            out = find_fwhm(imfile, star_pix[:,1][mask][mask2][mask3], star_pix[:,0][mask][mask2][mask3], plot=debug)
            mask_valid_fwhm = (out['detected']) * (out['e']>0.6) * ~np.isnan(out['fwhm']* (out['fwhm'] < 30))            

            if ((np.count_nonzero(mask_valid_fwhm) < 3) and (fitsutils.get_par(imfile, "FILTER")!="u")) or ( (np.count_nonzero(mask_valid_fwhm) < 2) and (fitsutils.get_par(imfile, "FILTER")=="u")):
                logger.error( "ERROR with FWHM!! Too few points for a valid estimation. %d"% np.count_nonzero(mask_valid_fwhm)+ ") points")
                logger.error( "%s %s"%(out["detected"], out["fwhm"]))
                return False

            outd = out[mask_valid_fwhm]

            logger.info( 'Average FWHM %.3f arcsec, %.3f pixels'%(np.median(outd['fwhm']),  np.median(outd['fwhm'])*pix2ang))
            
            fwhm = np.percentile(outd['fwhm'], 40)
            fitsutils.update_par(imfile,'FWHM', np.round(fwhm, 3))

            if (band in 'ugriz'):
                header='x y objid ra dec u g r i z du dg dr di dz'
            elif band in 'UBVRI':
                header='x y objid ra dec U B V R I dU dB dV dR dI'
            np.savetxt('/tmp/sdss_cat_det_%s.txt'%creationdate, z[mask_valid_fwhm], fmt=fmt, \
            header=header)
            print "Saved to", '/tmp/sdss_cat_det_%s.txt'%creationdate
            
            
        
    #Plot results
    img = img - np.nanmin(img)
    zmin, zmax = zscale.zscale(img)
        
    logger.info( "Found %d stars in %s. "%(len(cat_dec), survey)+\
	     "%d of them within the FoV. "%len(cat_ra[mask]) +\
            "%d of them are isolated."%len(cat_ra[mask][mask2])+\
            "%d of them with suitable magnitudes. "%len(cat_ra[mask][mask2][mask3]) +\
            "%d of them with detected stars."%np.count_nonzero(mask_valid_fwhm)) 
    
    
    if (plot):
        im = plt.imshow(img, aspect="equal", origin="lower", cmap=matplotlib.cm.gray_r, interpolation="none", vmin=zmin, vmax=zmax)

        
        if (len(star_pix[:,0][mask]) >0):        
            plt.scatter(star_pix[:,0][mask], star_pix[:,1][mask], marker="o", s=np.minimum(150, 10000*(10./mag[mask][mask2])**9), edgecolor="red", facecolor="none", label="catalogue")
        
        if (len(star_pix[:,0][mask][mask2]) >0):        
            plt.scatter(star_pix[:,0][mask][mask2], star_pix[:,1][mask][mask2], marker="o", s=20, edgecolor="yellow", facecolor="none", label="isolated")

        if (len(star_pix[:,0][mask][mask2][mask3]) >0):        
            plt.scatter(star_pix[:,0][mask][mask2][mask3], star_pix[:,1][mask][mask2][mask3], marker="o", s=200, edgecolor="green", facecolor="none", label="wihtin farme and mag")


        selected = star_pix[:,:][mask][mask2][mask3][mask_valid_fwhm]
        if (len(selected) >0):        
            plt.scatter(selected[:,0], selected[:,1], marker="o", \
                s=400, edgecolor="blue", facecolor="none")
            for i in np.arange(len(selected)):
                plt.text(selected[i,0]+10, selected[i,1]+10, i+1)
        
        plt.legend(loc="best", frameon=False, framealpha=0.9)
        plt.savefig( os.path.join( plotdir, os.path.basename(imfile).replace('.fits', '.seqstars.png').replace('.new', '.seqstars.png')))
        logger.info( "Saved stars to %s"%imfile.replace('.fits', '.seqstars.png'))
        plt.clf()
        
    return True
        
def add_to_zp_cal(ref_stars, image, logname):
    '''
    Records the parameters and instrumental and standard star magnitudes needed for zeropoint calibration.
    ref_stars: name of the file that contains the reference stars that are present in the image.
    image: fits image with the sources in them.
    logname: name of the file where the information on reference and measurements are logged.
    
    minst = M + c0 + c1(AIRMASS) + c2(color) + c3(UT)    
    '''
    
    coldic = {'u':'g', 'g':'r', 'r':'i', 'i':'r', 'z':'i', 'U':'B', 'B':'V', 'V':'R', 'R':'I', 'I':'R'}

    r = np.genfromtxt(ref_stars, delimiter=" ", dtype=None, names=True)
    imapp = os.path.join(os.path.join(os.path.dirname(image), "photometry"), os.path.basename(image) + ".app.mag")
    #imapp = os.path.join(os.path.dirname(image), os.path.basename(image) + ".app.mag")

    '''try:
        my = np.genfromtxt(imapp, comments="#", dtype=[("id","<f4"), ("X","<f4"), ("Y","<f4"),("Xshift","<f4"), ("Yshift","<f4"),("fwhm","<f4"), ("ph_mag","<f4"), ("stdev","<f4"), ("fit_mag","<f4"), ("fiterr","<f4")])
    except:'''
    my = np.genfromtxt(imapp, comments="#", dtype=[("id","<f4"), ("filename","<f4"),  ("X","<f4"), ("Y","<f4"),("Xshift","<f4"), ("Yshift","<f4"),("fwhm","<f4"), ("ph_mag","<f4"), ("stdev","<f4"), ("fit_mag","<f4"), ("fiterr","<f4")])

    if (my.size < 2):
        my = np.array([my])
    if (r.size < 2):
        r = np.array([r])
    
    band = fitsutils.get_par(image, 'filter')
    mask_valid1 = np.array(my["fwhm"]<9000) * np.array(my["ph_mag"]<9000) * np.array(~ np.isnan(r[band])) * np.array(~ np.isnan(my["fit_mag"]))
    
    r = r[mask_valid1]
    my = my[mask_valid1]
    N = len(r)
    
    my["fiterr"][np.isnan( my["fiterr"])] = 100

    col_band = coldic[band]
    exptime = fitsutils.get_par(image, "exptime")
    airmass = 1.3
    name = "object"
    if (fitsutils.has_par(image, "NAME")):
        name = fitsutils.get_par(image, "NAME")
    if (fitsutils.has_par(image, "AIRMASS")):
        airmass = fitsutils.get_par(image, "AIRMASS")
    if (fitsutils.has_par(image, "JD")):
        date = fitsutils.get_par(image, "JD")
    elif (fitsutils.has_par(image, "MJD")):
        date = fitsutils.get_par(image, "MJD")
    elif (fitsutils.has_par(image, "MJD-OBS")):
        date = fitsutils.get_par(image, "MJD-OBS")
    else:
        date = 0
    
    if (not os.path.isfile(logname)):
            with open( logname, "a") as f:
                f.write("#object,filename,filter,std,stderr,inst,insterr,jd,airmass,color,exptime\n")
    with open( logname, "a") as f:
        for i in range(N):
            f.write("%s,%s,%s,%.3f,%.3f,%3f,%.3f,%.3f,%.3f,%.3f,%.3f\n"%
            (name, image, band, r[i][band], r[i]["d"+band], my[i]['fit_mag'], my[i]['fiterr'], date, airmass, r[i][band]-r[i][col_band],exptime))

def lsq_test():

    x = np.random.normal(0, 1, 100)
    y = np.random.normal(0, 1, 100)
    
    X, Y = np.meshgrid(x, y)
    
    D = 23 + 2*X + 4*Y + np.random.rand(100, 100)*1
    O = np.zeros_like(D)
    
    M = np.zeros((len(X.flatten()), 3))
    M[:,0] = 1 
    M[:,1] = X.flatten()
    M[:,2] = Y.flatten()
    
    depend = D.flatten() - O.flatten()
        
    lsq_result = np.linalg.lstsq(M, depend)
    coef = lsq_result[0]
    res = lsq_result[1]
    
    print coef, res
    
    f, (ax1, ax2) = plt.subplots(2, 1)
    ax1.plot(X.flatten(), depend -coef[0] -Y.flatten()*coef[2] , ".")
    ax1.plot(X.flatten(),  X.flatten()*coef[1])
    ax1.set_xlabel("X")

    ax2.plot(Y.flatten(), depend -coef[0] -X.flatten()*coef[1] , ".")
    ax2.plot(Y.flatten(),  Y.flatten()*coef[2])
    ax2.set_xlabel("Y")
    
    plt.show()
       
       
def calibrate_zp_fourshot(logfile, plot=True):
    '''
    The field of view is quite small, therefore, all rc shots are used to calibrate the zeropoint.
    This routine retrieves all the stars taken with the same filter when pointing to the science object.
    '''

    a = np.genfromtxt(logfile, dtype=None, names=True, delimiter=",")
    a.sort(order=['jd'], axis=0)
    a = a[a['inst']!=0]

    plotdir = os.path.join(os.path.dirname(os.path.abspath(logfile)), "photometry")

    for name in set(a['object']):    
        for b in set(a['filter']):
            aib = a[(a["filter"]==b)*(a["object"]==name)*(np.abs(a["color"])<1)]
            
            if len(aib) < 3:
                print "Less than 3 stars found with %s %s"%(name,b)
                continue
            
            #First fit for the linear and detect the deviations
            coefs, residuals, rank, singular_values, rcond = np.polyfit(aib["std"], aib["inst"], w=1./np.maximum(0.3, np.sqrt(aib["stderr"]**2 + aib["insterr"]**2)), deg=1, full=True)
            p = np.poly1d(coefs)
            
            if (plot):
                plt.figure()
                plt.title("%s Filter %s"%(name, b))
                plt.errorbar(aib["std"], aib["inst"], yerr=np.sqrt(aib["stderr"]**2 + aib["insterr"]**2), fmt="o")
                plt.plot(aib["std"], p(aib["std"]))
                
            diff = np.abs(aib["inst"] - p(aib["std"]))
            mad = stats.funcs.median_absolute_deviation(diff)
            aib = aib[diff<mad*5]
            
            if (plot):
                plt.errorbar(aib["std"], aib["inst"], yerr=np.sqrt(aib["stderr"]**2 + aib["insterr"]**2), fmt="o")
                plt.plot(aib["std"], p(aib["std"]))
                plt.savefig(os.path.join(plotdir, "zp_mag_mag_%s_%s.png"%(name, b)))
                plt.close()
            

            #Then fit for the colour
            coefs, residuals, rank, singular_values, rcond = np.polyfit(aib["color"], aib["std"] - aib["inst"], w=1./np.maximum(0.15, np.sqrt(aib["stderr"]**2 + aib["insterr"]**2)), deg=1, full=True)
            p = np.poly1d(coefs)
            
            color, zp = coefs
            
            mad = stats.funcs.median_absolute_deviation(p(aib["color"]) - (aib["std"] - aib["inst"]))
             
            print coefs, residuals, mad
            
            for f in aib["filename"]:
                #Add these values to the header.
                pardic = {"IQZEROPT" : 1,\
                        "ZPCAT" : "SDSS4shot",\
                        "ZEROPTU" : float("%.3f"%mad),\
                        "ZEROPT" : float("%.3f"%zp),\
                        "ZP": float("%.3f"%zp),\
                        "ZPERR": float("%.3f"%mad)}
                fitsutils.update_pars(f, pardic)

            if (plot):
                plt.figure()
                plt.title("%s Filter %s"%(name, b))
                plt.errorbar(aib["color"], aib["std"] - aib["inst"], yerr=np.sqrt(aib["stderr"]**2 + aib["insterr"]**2), fmt="o")
                plt.plot(aib["color"], p(aib["color"]))
                plt.savefig(os.path.join(plotdir, "zp_col_mag_%s_%s.png"%(name, b)))
                plt.close()
    return
    
def lsq_zeropoint(logfile, plotdir=None, plot=True):
    '''
    Uses least squares approach to compute the optimum coefficients for ZP, colour term, airmass and time.
    
    '''
    a = np.genfromtxt(logfile, dtype=None, names=True, delimiter=",")
    a.sort(order=['jd'], axis=0)
    a = a[a['inst']!=0]
    
    a['jd'] = a['jd'] - np.min(a['jd'])
    cols = {'u':'purple', 'g':'green', 'r':'red', 'i':'orange'}



    for b in set(a['filter']):
        ab = a[a['filter']==b]
        
        #Filter extreme colours which may bias the coefficient.
        mincol = np.median(ab["color"]) - 2*np.std(ab["color"])
        maxcol = np.median(ab["color"]) + 2*np.std(ab["color"])
        
        ab = ab[ (ab["color"]>mincol) * (ab["color"]<maxcol) ]        
        
        #Remove detections which are too far
        coefs, residuals, rank, singular_values, rcond = np.polyfit(ab["std"], ab["inst"], w=1./np.maximum(0.3, np.sqrt(ab["stderr"]**2 + ab["insterr"]**2)), deg=1, full=True)
        p = np.poly1d(coefs)
        
        if (plot):
            plt.figure()
            plt.title("Filter %s"%(b))
            plt.errorbar(ab["std"], ab["inst"], yerr=np.sqrt(ab["stderr"]**2 + ab["insterr"]**2), fmt="o")
            plt.plot(ab["std"], p(ab["std"]))
            
        diff = np.abs(ab["inst"] - p(ab["std"]))
        mad = stats.funcs.median_absolute_deviation(diff)
        ab = ab[diff<mad*5]
        
                
        if (plot):
            plt.figure()
            plt.title("Filter %s"%(b))
            plt.errorbar(ab["std"], ab["inst"], yerr=np.sqrt(ab["stderr"]**2 + ab["insterr"]**2), fmt="o")
            plt.plot(ab["std"], p(ab["std"]))
            plt.show()
            
        #Find the coefficients.
        '''M = np.zeros((len(ab), 4))
        M[:,0] = 1      
        #M[:,1] = ab['inst']
        M[:,1] = ab['color'] 
        M[:,2] = ab['airmass'] - 1.3
        M[:,3] = ab['jd'] 
        
        #print M
        
        depend = ab['std']-ab['inst']
        
        lsq_result = np.linalg.lstsq(M, depend)
        coef = lsq_result[0]
        res = lsq_result[1]
                
        print "Band", b, "zp %.2f, col %.2f , airmass %.3f , time %.2f"%(coef[0], coef[1], coef[2], coef[3]), "Residuals", res
        
        emp_col = depend -coef[0] -(ab['airmass']-1.3)*coef[2] - ab['jd']*coef[3]
        pred_col = ab['color']*coef[1]
        mask_col_outlier = np.abs(emp_col - pred_col) > 3* stats.funcs.median_absolute_deviation(emp_col - pred_col)


        emp_airmass = depend -coef[0] - ab['color']*coef[1] - ab['jd']*coef[3]
        pred_airmass = ab['airmass']*coef[1]
        mask_airmass_outlier = np.abs(emp_airmass - pred_airmass) > 3 * stats.funcs.median_absolute_deviation(emp_airmass - pred_airmass)
        
        emp_jd = depend -coef[0] -ab['color']*coef[1]- (ab['airmass']-1.3)*coef[2]
        pred_jd =  ab['jd']*coef[3]
        mask_col_jd = np.abs(emp_jd - pred_jd) > 3 * stats.funcs.median_absolute_deviation(emp_jd - pred_jd)
        
        mask = mask_col_outlier * mask_airmass_outlier * mask_col_jd
        
        ab = ab[~mask]'''
        M = np.zeros((len(ab), 5))
        M[:,0] = 1      
        #M[:,1] = ab['inst']
        M[:,1] = ab['color'] 
        M[:,2] = ab['airmass'] - 1.3
        M[:,3] = ab['jd'] 
        M[:,4] = ab['jd']**2
        #M[:,5] = ab['jd']**3
        
        #print M
        
        depend = ab['std']-ab['inst']
        
        lsq_result = np.linalg.lstsq(M, depend)
        coef = lsq_result[0]
        res = lsq_result[1]
        
        #Save the coefficients 
        np.savetxt("coefs_%s.txt"%b, coef)
        
        #Empirical and predicted values
        emp_col = depend -coef[0] -(ab['airmass']-1.3)*coef[2] - (coef[3]*ab['jd'] + coef[4]*ab['jd']**2)# + coef[5]*ab['jd']**3 )
        pred_col = ab['color']*coef[1]

        emp_airmass = depend -coef[0] - ab['color']*coef[1] - ( coef[3]*ab['jd'] + coef[4]*ab['jd']**2)# + coef[5]*ab['jd']**3)
        pred_airmass = (ab['airmass']-1.3)*coef[2]
        
        emp_jd = depend -coef[0] -ab['color']*coef[1]- (ab['airmass']-1.3)*coef[2]
        pred_jd =  coef[3]*ab['jd'] + coef[4]*ab['jd']**2# + coef[5]*ab['jd']**3 
                
        est_zp = coef[0] +ab['color']*coef[1] +(ab['airmass']-1.3)*coef[2] + coef[3]*ab['jd'] + coef[4]*ab['jd']**2 #+ coef[5]*ab['jd']**3 
        rms =  np.sqrt(np.sum((depend-est_zp)**2)/(len(depend)-1))
        np.savetxt("rms_%s.txt"%b, np.array([rms]))
        print "Filter %s ZP %.2f RMS %.2f"%(b, coef[0], rms)

        if (plot):
            plt.close("all")
            f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
            ax1.plot(ab['color'], emp_col, "o", color=cols[b], ms=4, alpha=0.4)
            ax1.plot(ab['color'],  pred_col, color=cols[b])
            ax1.set_xlabel("color")

            ax2.plot(ab['airmass'], emp_airmass, "o", color=cols[b], ms=4, alpha=0.4)
            ax2.plot(ab['airmass'], pred_airmass , color=cols[b])
            ax2.set_xlabel("airmass")

            #arr = np.array([ab['jd'], pred_jd]).T
            #print arr, arr.shape, arr.dtype
            #arr.view(dtype=[('f0', np.float64), ('f1', np.float64)]).sort(order=['f0'], axis=0)
            
            ax3.plot(ab['jd'], emp_jd, "o", color=cols[b], ms=4, alpha=0.4)
            ax3.plot(ab['jd'], pred_jd, color=cols[b])
            ax3.set_xlabel("jd")   
            
            #ax4.plot(ab['jd'], depend, "*", color=cols[b], alpha=0.4)
            #ax4.errorbar(ab['jd'], est_zp, yerr=ab['insterr'], marker="o", c=cols[b], ls="none")
            ax4.errorbar(ab['jd'], depend-est_zp, yerr=np.minimum(1, ab['insterr']), fmt="o", c=cols[b], alpha=0.4, ms=4)
            ax4.set_xlabel("jd")  
            ax4.set_ylabel("night predicted zeropoint")
            ax4.invert_yaxis()
            
            if (not plotdir is None):
                plt.savefig(os.path.join(plotdir, "allstars_%s.png"%b))
            else:	
                plt.show()

def interpolate_zp(reduced, logfile):
    '''
    Uses the zeropoint coefficients derived from SDSS fields to interpolate the 
    zeropoint for the images that are outside of SDSS field.
    '''        
    a = np.genfromtxt(logfile, dtype=None, names=True, delimiter=",")
    a.sort(order=['jd'], axis=0)
    a = a[a['inst']!=0]
    a = a[a['insterr']>0]
    
    jdmin =  np.min(a['jd'])
    
    jdmax = np.max(a['jd']) - jdmin
    
    zpfiles = glob.glob(os.path.join(reduced, "*fits"))
    
    zpfiles = [zf for zf in zpfiles if fitsutils.has_par(zf, "IQZEROPT") and \
    (fitsutils.get_par(zf, "IQZEROPT")==0 or fitsutils.get_par(zf, "ZEROPT")==0)]
    
    #Load the coefficients.
    coefs = {}    
    for fi in ["u", "g", "r", "i"]:
        coeffile = os.path.join(reduced, "coefs_%s.txt"%fi)
        if os.path.isfile(coeffile):
            coefs[fi] = np.genfromtxt(coeffile)

    #Load the rms.
    rms = {}    
    for fi in ["u", "g", "r", "i"]:
        rmsfile = os.path.join(reduced, "rms_%s.txt"%fi)
        if os.path.isfile(rmsfile):
            rms[fi] = np.genfromtxt(os.path.join(reduced, "rms_%s.txt"%fi))
        
    for image in zpfiles:
        filt = fitsutils.get_par(image, "FILTER")
        
        #To not extrapolate outside of the valid interval.
        jd = np.minimum(jdmax, fitsutils.get_par(image, "JD") - jdmin)
        airmass = fitsutils.get_par(image, "AIRMASS")
        
        #If there are coefficients for that filter, load them and interpolate.
        #Otherwise, skip this file.
        if not coefs.has_key(filt):
            continue
        
        #est_zp = coef[0] +ab['color']*coef[1] +(ab['airmass']-1.3)*coef[2] + coef[3]*ab['jd'] + coef[4]*ab['jd']**2 + coef[5]*ab['jd']**3 + coef[6]*ab['jd']**4 + coef[7]*ab['jd']**5

        values = np.array([1, 0, airmass-1.3, jd, jd**2])
        est_zp = np.sum(coefs[filt]*values)
        
        #Update the header with the computed zeropoint.
        pardic = {
                "IQZEROPT" : 1,\
                "ZPCAT" : "SDSSinterpolated",\
                "ZEROPTU" : float(rms[filt]),\
                "ZEROPT" : est_zp}
        fitsutils.update_pars(image, pardic)
    
    
def lsq_zeropoint_partial(logfile, plot=True):
    '''
    Uses least squares approach to compute the optimum coefficients for ZP, colour term, airmass and time.
    
    '''
    a = np.genfromtxt(logfile, dtype=None, names=True, delimiter=",")
    
    a = a[a['insterr']<0.1]
    a['jd'] = a['jd'] - np.min(a['jd'])
    cols = {'u':'purple', 'g':'green', 'r':'red', 'i':'orange'}


    for b in set(a['filter']):
        ab = a[a['filter']==b]
        M = np.zeros((len(ab), 2))
        M[:,0] = 1      
        #M[:,1] = ab['inst']
        M[:,1] = ab['color'] 

        
        #print M
        
        depend = ab['std']-ab['inst']
        
        lsq_result = np.linalg.lstsq(M, depend)
        coef = lsq_result[0]
        res = lsq_result[1]
                
        print "Band", b, "zp %.2f, col %.2f"%(coef[0], coef[1]), "Residuals", res
        
        if (plot):
            f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
            ax1.plot(ab['color'], depend -coef[0], "o", color=cols[b])
            ax1.plot(ab['color'],  ab['color']*coef[1] , color=cols[b])
            ax1.set_xlabel("color")
     
    if (plot):
        plt.show()
        
def find_zeropoint_noid(ref_stars, image, plot=True, plotdir="."):
    '''
    Finds the zeropoint by comparig the magnitude of the stars measured in the image,
    vs. the magnitude of the stars from the reference catalogue.
    band: band for which we want to measure the zeropoint.
    col_band: band for the colour term we want to use.
    
    returns the zeropoint, the colour term and the standard deviation in the 
    zeropoint from all the measured stars.
    
    '''
    
    logger.info( 'Finding the optimum ZP fit...')
    
    logger.info("Reference stars used: %s"% ref_stars)
    
    r = np.genfromtxt(ref_stars, delimiter=" ", dtype=None, names=True)
    imapp = os.path.join(os.path.join(os.path.dirname(image), "photometry"), os.path.basename(image) + ".app.mag")
    #imapp = os.path.join(os.path.dirname(image), os.path.basename(image) + ".app.mag")

    my = np.genfromtxt(imapp, comments="#", dtype=[("id","<f4"), ("image","|S20"),   ("X","<f4"), ("Y","<f4"),("Xshift","<f4"), ("Yshift","<f4"),("fwhm","<f4"), ("ph_mag","<f4"), ("stdev","<f4"), ("fit_mag","<f4"), ("fiterr","<f4")])
    '''try:
    except:
        my = np.genfromtxt(imapp, comments="#", dtype=[("id","<f4"), ("X","<f4"), ("Y","<f4"),("Xshift","<f4"), ("Yshift","<f4"),("fwhm","<f4"), ("ph_mag","<f4"), ("stdev","<f4"), ("fit_mag","<f4"), ("fiterr","<f4")])
    '''
    
    if (my.size < 2):
        my = np.array([my])
    if (r.size < 2):
        r = np.array([r])
    coldic = {'u':'g', 'g':'r', 'r':'i', 'i':'z', 'z':'i', 'U':'B', 'B':'V', 'V':'R', 'R':'I', 'I':'R'}
    band = fitsutils.get_par(image, 'filter')
    col_band = coldic[band]
    
    mask_valid1 = np.array(my["fwhm"]<9000) * np.array(my["ph_mag"]<9000) * np.array(~ np.isnan(r[band])) * np.array(~ np.isnan(my["fit_mag"]))
    
    N = len(r)
    r = r[mask_valid1]
    my = my[mask_valid1]
    
    if len(my) == 0:
        logger.warn( "Warning, no reliable measurements for file %s. Returning 0,0,0."% image)
        return 0, 0, 0
        
    my["fiterr"][np.isnan( my["fiterr"])] = 100
    
    ids = np.arange(N)+1
    ids = ids[mask_valid1]
    
    coefs, residuals, rank, singular_values, rcond = np.polyfit(r[band]-r[col_band], r[band] - my["fit_mag"], w=1./np.maximum(0.1, np.sqrt(my["fiterr"]**2 + r['d'+band]**2)), deg=1, full=True)
    p = np.poly1d(coefs)
    
    logger.info("Coefficients for 1 deg polynomial fit to the zeropoint: %s"% coefs)
    logger.info("%s - %s = %.3f, %.3f"%(band, col_band, p[0], p[1]))


    pred = p(r[band]-r[col_band])
    measured = r[band]- my["fit_mag"]
    mad = stats.funcs.median_absolute_deviation(pred-measured)
    
    print "MAD1",  mad
    
    if len(r) > 4:
        mask = np.abs(pred-measured)/mad < 3
            
        r = r[mask]
        my = my[mask]
        ids = ids[mask]
        
        coefs, residuals, rank, singular_values, rcond = np.polyfit(r[band]-r[col_band], r[band] - my["fit_mag"], w=1./np.maximum(0.2, np.sqrt(my["fiterr"]**2 + r['d'+band]**2)), deg=1, full=True)
        
        logger.info("Coefficients for 1 deg polynomial fit to the zeropoint: %s. After outlier rejection."% coefs)
        
        p = np.poly1d(coefs)
    
        pred = p(r[band]-r[col_band])
        measured = r[band]- my["fit_mag"]
        mad = stats.funcs.median_absolute_deviation(pred-measured)
        print "MAD2",  mad

        
    if (plot):
        print "Plotting..."
        plt.errorbar(r[band]-r[col_band], r[band] - my["fit_mag"] , yerr=np.sqrt(my["fiterr"]**2 + r['d'+band]**2), marker="o", ls="None")
        for i, myid in enumerate(ids):
            plt.text(r[band][i]-r[col_band][i] + 0.01, r[band][i] - my["fit_mag"][i]+ 0.01, str(myid))
        x = np.linspace(np.min(r[band]-r[col_band]), np.max(r[band]-r[col_band]), 100)
        plt.plot(x, p(x))
        plt.title("Best fit ZP: %.2f colour term: %.2f MAD: %.2f"%(p[0], p[1], mad))
        plt.xlabel("{:} - {:}".format(band, col_band))
        plt.ylabel("ZP")
        plt.savefig( os.path.join( plotdir, os.path.basename(image).replace('.fits', ".zp.png").replace('.new', ".zp.png")))

        plt.clf()
    
    logger.info("%s - %s = %.3f, %.3f"%(band, col_band, p[0], p[1]))
    
    pred = p(r[band]-r[col_band])
    measured = r[band]- my["fit_mag"]

    fitsutils.update_par(image, "ZP", np.round(p[0], 3))
    fitsutils.update_par(image, "COLTERM", np.round(p[1], 3))
    fitsutils.update_par(image, "ZPERR", np.round(mad, 3)) #np.std(pred - measured))


    return np.round(p[0], 3), np.round(p[1], 3), np.round(mad, 3)


def calibrate_zeropoint(image, plot=True, plotdir=None, debug=False, refstars=None):
    '''
    Calibrates the zeropoint using SDSS catalogue.    
    '''
    
    
    if plot and plotdir is None:
        plotdir = os.path.join(os.path.dirname(image), "photometry")
        if not os.path.isdir(plotdir):
            os.makedirs(plotdir)
        
    filt = fitsutils.get_par(image, 'filter')
    exptime = fitsutils.get_par(image, "exptime")
    if fitsutils.has_par(image, "JD"):
        date = fitsutils.get_par(image, "JD")
    elif fitsutils.has_par(image, "MJD"):
        date = fitsutils.get_par(image, "MJD")
    elif fitsutils.has_par(image, "MJD-OBS"):
        date = fitsutils.get_par(image, "MJD-OBS")
    else:
        date=0
        
    if fitsutils.has_par(image, "AIRMASS"): 
        airmass = fitsutils.get_par(image, "AIRMASS")
    else:
        airmass = 1.3
    objname = fitsutils.get_par(image, "OBJECT")
    band = fitsutils.get_par(image, "FILTER")
    
    logger.info( "Starting calibration of ZP for image %s for object %s with filter %s."%(image, objname, band))

    if (exptime < 10):
        logger.error( "ERROR. Exposure time too short for image (%s) to see anything..."%image)
        return 

    extracted = extract_star_sequence(os.path.abspath(image), filt, plot=plot, survey='sdss', debug=debug, refstars=refstars, plotdir=plotdir)
    if (not extracted):
        logger.warn( "Field not in SDSS or error when retrieving the catalogue... Skipping. Image %s not zeropoint calibrated."%image)
            #Add these values to the header.
        pardic = {"IQZEROPT" : 0,\
            "ZPCAT" : "None",\
            "ZEROPTU" : 0,\
            "ZEROPT" : 0, \
            "ZP":0,\
            "ZPERR":0}
        fitsutils.update_pars(image, pardic)
        return 

    #If extraction worked, we can get the FWHM        
    fwhm = fitsutils.get_par(image, "fwhm")
    fwhm_as = fwhm * 0.394

    app_phot.get_app_phot("/tmp/sdss_cat_det_%s.txt"%creationdate, image, wcsin='logic', plotdir=plotdir, box=20)
    
    #Compute the zeropoint for the specific image.
    z, c, err = find_zeropoint_noid("/tmp/sdss_cat_det_%s.txt"%creationdate, image, plot=plot, plotdir=plotdir)
    
    #Add these values to the header.
    pardic = {"IQZEROPT" : 1,\
            "ZPCAT" : "SDSS",\
            "ZEROPTU" : np.round(err, 3),\
            "ZEROPT" : np.round(z, 3)}
    fitsutils.update_pars(image, pardic)
            
    #Log the current zeropoint for this image
    logname = os.path.join(os.path.dirname(image), "zeropoint.log")
    
    #Add the data to a later stage zeropoint calibrtion with all-sky data.
    zplogname = os.path.join(os.path.dirname(image), "allstars_zp.log")
    
    add_to_zp_cal("/tmp/sdss_cat_det_%s.txt"%creationdate, image, zplogname)


    if (not os.path.isfile(logname)):
            with open( logname, "a") as f:
                f.write("#filename,exptime,filter,date,airmass,fwhm_pix,fwhm_as,zeropoint,color,err\n")
    with open( logname, "a") as f:
        f.write("%s,%.1f,%s,%3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n"%(image,exptime,filt,date,airmass,fwhm,fwhm_as,z,c,err))
        
    clean_tmp_files()
        
def plot_zp(zpfile, plotdir=None):
    import datetime
    
    def mkdate(text):
        return datetime.datetime.strptime(text, '%Y-%m-%dT%H:%M:%S') 
    

    plt.clf()
    a = np.genfromtxt(zpfile, names=True, dtype=None, delimiter=',')
    
    cols = {'u':'purple', 'g':'green', 'r':'red', 'i':'orange'}
    
    for fi in ['u', 'g', 'r', 'i']:
        for i in range(len(a[a['filter']==fi])):
            plt.errorbar( (a[a['filter']==fi]['date'][i]-np.min(a[a['filter']==fi]['date'], axis=0))*24., \
            a[a['filter']==fi]['zeropoint'][i], yerr=a[a['filter']==fi]['err'][i], marker='o', mfc=cols[fi], mec='k', ecolor=cols[fi], ls='none', ms=a[a['filter']==fi]['fwhm_pix'][i])
        logger.info( "Median zeropoint for filter %s: %.2f mag"%(fi, np.median(a[a['filter']==fi]['zeropoint'])))

    plt.gca().invert_yaxis()
    plt.xlabel("Obs Date (JD - min(JD)) [h]")
    plt.ylabel("ZP [mag]")
    plt.ylim(25,18)
    if (plotdir is None):
        plt.show()
    else:
	plt.savefig(os.path.join(plotdir, "zeropoint_per_exposure.png"))
	plt.clf()
    
def plot_zp_airmass(zpfile):
    import datetime
    
    def mkdate(text):
        return datetime.datetime.strptime(text, '%Y-%m-%dT%H:%M:%S') 
    
    a = np.genfromtxt(zpfile, names=True, dtype=None, delimiter=',')
    
    cols = {'u':'purple', 'g':'green', 'r':'red', 'i':'orange'}
    
    for fi in ['u', 'g', 'r', 'i']:
        for i in range(len(a[a['filter']==fi])):
            plt.errorbar( a[a['filter']==fi]['airmass'][i], a[a['filter']==fi]['zeropoint'][i], yerr=a[a['filter']==fi]['err'][i], marker='o', mfc=cols[fi], mec='k', ecolor=cols[fi], ls='none', ms=a[a['filter']==fi]['fwhm_pix'][i])
        logger.info( "Median zeropoint for filter %s: %.2f mag"%(fi, np.median(a[a['filter']==fi]['zeropoint'])))

    plt.gca().invert_yaxis()
    plt.xlabel("Airmass")
    plt.ylabel("ZP [mag]")
    plt.show()


def reset_zp(directory):
    '''
    Resets the zeropoint keyword to if the result is interpolated.
    '''    
    
    files = glob.glob(directory+"/*fits")
    
    for f in files:
        if (fitsutils.get_par(f, "ZPCAT")=="SDSSinterpolated"):
            fitsutils.update_par(f, "IQZEROPT", 0)
            fitsutils.update_par(f, "ZEROPT", 0)
            fitsutils.update_par(f, "ZEROPTU", 0)
    
    
def main(reduced):
    '''
    Performs the main zeropoint calculations for the folder and plots the results.
    
    '''
    os.chdir(reduced)
    
    plotdir = "zeropoint"
    if (not os.path.isdir(plotdir)):
        os.makedirs(plotdir)
    

    for f in glob.glob("*.fits|*.new"):
        logger.info("Starting calibration of zeropoint for %s"% f)
        if (not fitsutils.has_par(f, "IMGTYPE") or fitsutils.get_par(f, "IMGTYPE") == "SCIENCE"):
            calibrate_zeropoint(f, plotdir=os.path.abspath(plotdir))
    if (os.path.isfile("zeropoint.log")):
        plot_zp("zeropoint.log", plotdir)
    if (os.path.isfile("allstars_zp.log")):
        lsq_zeropoint("allstars_zp.log", plotdir)
        interpolate_zp(reduced, "allstars_zp.log")
        
    clean_tmp_files()
     
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description=\
        '''

        Runs astrometry.net on the image specified as a parameter and returns 
        the offset needed to be applied in order to center the object coordinates 
        in the reference pixel.
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('reduced', type=str, help='Directory containing the reduced fits for the night.')

    args = parser.parse_args()
    
    reduced = args.reduced
    
    if (reduced is None):
        print "Please, add the directory containing reduced data as a parameter."
    else:
        main(reduced)

