# -*- coding: utf-8 -*-
"""
Created on Thu Dec  3 14:57:05 2015

@author: nadiablago
"""
import zscale
import sys, urllib
import pyfits as pf
import numpy as np
import matplotlib.pyplot as plt
import math
import pywcs 
from astropy.io.votable import parse_single_table
import coordinates_conversor
import matplotlib
import app_phot
from numpy.lib import recfunctions as rfn
import scipy.optimize as opt
import fitsutils
import os, shutil

def are_isolated(rav, decv, r):
    '''
    Returns a mask saying whether the stars are isolated in a radius of r arcsec.
    '''
    
    index = np.arange(len(rav))
    mask = []
    
    for i in np.arange(len(rav)):
        d = coordinates_conversor.get_distance(rav[i], decv[i], rav[index!=i], decv[index!=i])
        mask.append(~ np.any(d*3600<r))
    
    return np.array(mask)

def twoD_Gaussian(xdata_tuple, amplitude, xo, yo, sigma_x, sigma_y, theta, offset):
    (x, y) = xdata_tuple                                                        
    xo = float(xo)                                                              
    yo = float(yo)                                                              
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)   
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)    
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)   
    g = offset + amplitude*np.exp( - (a*((x-xo)**2) + 2*b*(x-xo)*(y-yo) + c*((y-yo)**2)))                                   
    return g.ravel()
    
def twoD_Gauss_test():
    # Create x and y indices
    x = np.linspace(0, 200, 201)
    y = np.linspace(0, 200, 201)
    x, y = np.meshgrid(x, y)
    
    #create data
    data = twoD_Gaussian((x, y), 3, 100, 100, 20, 40, 0, 10)
    
    # plot twoD_Gaussian data generated above
    plt.figure()
    plt.imshow(data.reshape(201, 201))
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
        hrad = math.ceil(rad/2)
        sub = img[x_i-hrad:x_i+hrad, y_i-hrad:y_i+hrad]
        x = np.linspace(0, len(sub), len(sub))
        y = np.linspace(0, len(sub), len(sub))
        x, y = np.meshgrid(x, y)
        print x.shape, sub.shape
    
        #(xdata_tuple, amplitude, xo, yo, def_fwhm, def_fwhm, theta, offset):
        def_x = np.argmax(np.sum(sub, axis=0))
        def_y = np.argmax(np.sum(sub, axis=1))

        initial_guess = (100, def_x, def_y, def_fwhm, def_fwhm, 0, np.percentile(sub, 50))
        detected = True
        try:
            popt, pcov = opt.curve_fit(twoD_Gaussian, (x, y), sub.flatten(), p0=initial_guess)
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
        
        
        print i, detected, 'Amplitude %.3f\t'%amplitude, 'BG %.3f\t'%background, "BG_stats %.3f\t"%np.percentile(sub, 50), "FWHM_x,FWHM_y=(%.3f, %.3f)"%(fwhm_x, fwhm_y)
        
        out[i] = (detected, np.average([fwhm_x, fwhm_y]), np.minimum(fwhm_x, fwhm_y) / np.maximum(fwhm_x, fwhm_y))
        
        if (detected & plot):
            data_fitted = twoD_Gaussian((x, y), *popt)
            
            fig, (ax, ax2) = plt.subplots(1, 2)
            ax.hold(True)
            ax.imshow(sub, cmap=plt.cm.jet, origin='bottom', extent=(x.min(), x.max(), y.min(), y.max()))
            ax.contour(x, y, data_fitted.reshape(sub.shape[0], sub.shape[1]), 5, colors='w')
            ax2.imshow(sub-data_fitted.reshape(sub.shape[0], sub.shape[1]), cmap=plt.cm.jet, origin='bottom', extent=(x.min(), x.max(), y.min(), y.max()))
            ax2.contour(x, y, data_fitted.reshape(sub.shape[0], sub.shape[1]), 5, colors='w')
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
    
    
def extract_star_sequence(imfile, band, plot=True, survey='apass', debug=False, refstars=None):
    '''
    Plots the USNOB1 field of stars on top of the star field image.
    Red circles are stars identified from the catalogue in the correct magnitude range.
    Yellow circles are stars that are isolated.
    
    '''
    
    survey = str.lower(survey)
    
    f = pf.open(imfile)
    wcs = pywcs.WCS(f[0].header)

    img = f[0].data
    img[img<0] = 0
    
    ra, dec = wcs.wcs_pix2sky(np.array([img.shape[0]/2, img.shape[1]/2], ndmin=2), 1)[0]
    sr = 5.5/60
    print ra,dec
    
    if not refstars is None:
        shutil.copy(refstars, "/tmp/tmp_sdss.cat")
        catalog = np.genfromtxt("/tmp/tmp_sdss.cat", names=True, dtype=None, delimiter=",")
        cat_ra = catalog["ra"]
        cat_dec = catalog["dec"]
        mag = catalog["R"]
        
    elif str.lower(survey) =='usnob1':
        #ra, dec = coordinates_conversor.hour2deg(f[0].header['RA'], f[0].header['DEC'])
        #SEDM FoV is 6.5 arcmin, due to uncertainties in the position, 4 arcmin radius assumed.
        # Download USNO-B1 catalog for the position
        catalog_url = 'http://www.nofs.navy.mil/cgi-bin/vo_cone.cgi?CAT=USNO-B1&RA=%.5f&DEC=%.5f&SR=%.4f&VERB=1' % (ra, dec, sr)
        print "Downloading USNO-B1 catalog..."
        urllib.urlretrieve(catalog_url, '/tmp/tmp.cat')
        
        # Read RA, Dec and magnitude from XML format USNO catalog
        catalog = parse_single_table("/tmp/tmp.cat")
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
        urllib.urlretrieve(catalog_url, '/tmp/tmp_apass.cat')
        catalog = np.genfromtxt("/tmp/tmp_apass.cat", delimiter=",", names=True)
        
        cat_ra = catalog['radeg']
        cat_dec = catalog['decdeg']
        mag = catalog['Sloan_r']
        
    elif (survey=='sdss'):
        minmag = 0
        maxmag = 21.0
        catalog_url='http://skyserver.sdss.org/dr7/en/tools/search/x_radial.asp?ra=%.5f&dec=%.5f&check_type=type&type=6\
        &radius=%.4f&check_u=u&min_u=%.2f&max_u=%.2f&check_g=g&min_g=%.2f&max_g=%.2f&check_r=r&min_r=%.2f&max_r=%.2f&check_i=i&min_i=%.2f&max_i=%.2f&check_z=z&min_z=%.2f&max_z=%.2f&entries=top&topnum=500&format=csv'%(ra, dec, sr*60,minmag,maxmag,minmag,maxmag,minmag,maxmag,minmag,maxmag,minmag,maxmag)
        print "Downloading SDSS catalog..."
        print catalog_url
        urllib.urlretrieve(catalog_url, '/tmp/tmp_sdss.cat')
        catalog = np.genfromtxt("/tmp/tmp_sdss.cat", delimiter=",", names=True)

        try:
            cat_ra = catalog['ra']
            cat_dec = catalog['dec']
            mag = catalog[band]
        except:
            print "Problems with SDSS image", band
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
    mask2 = np.array(are_isolated(cat_ra[mask], cat_dec[mask], 30.))
        
    #Select only stars that are within the proper magnitude range
    mask3 = (mag[mask][mask2] < 20.) * (mag[mask][mask2] > 12) 
    
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
        if not np.any(mask) and not np.any(mask2) and not np.any(mask3):
            print star_pix
            print "No stars left...", mask, mask2, mask3
        else:
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
    
            np.savetxt('/tmp/sdss_cat.txt', z, fmt=fmt, \
            header='x y objid run rerun camcol field obj type ra dec u g r i z du dg dr di dz')
            print "Saved to", '/tmp/sdss_cat.txt'
            
            #Find FWHM for this image            
            out = find_fwhm(imfile, star_pix[:,1][mask][mask2][mask3], star_pix[:,0][mask][mask2][mask3], plot=debug)
            mask_valid_fwhm = (out['detected']) * (out['e']>0.7) * ~np.isnan(out['fwhm']* (out['fwhm'] < 30))            

            if (np.count_nonzero(mask_valid_fwhm) < 3):
                print "ERROR with FWHM!! Too few points for a valid estimation."
                return False

            outd = out[mask_valid_fwhm]

            print 'Average FWHM',np.median(outd['fwhm']), 'arcsec', np.median(outd['fwhm'])*pix2ang, 'pixels'
            
            fwhm = np.median(outd['fwhm'])
            fitsutils.update_par(imfile,'FWHM',fwhm)


            np.savetxt('/tmp/sdss_cat_det.txt', z[mask_valid_fwhm], fmt=fmt, \
            header='x y objid run rerun camcol field obj type ra dec u g r i z du dg dr di dz')
            print "Saved to", '/tmp/sdss_cat_det.txt'
            
            
        
    #Plot results
    img = img - np.nanmin(img)
    zmin, zmax = zscale.zscale(img)
        
    print "Found %d stars in %s. "%(len(cat_dec), survey), \
            "%d of them within the FoV. "%len(cat_ra[mask]),\
            "%d of them are isolated."%len(cat_ra[mask][mask2]),\
            "%d of them with suitable magnitudes. "%len(cat_ra[mask][mask2][mask3]),\
            "%d of them with detected stars."%np.count_nonzero(mask_valid_fwhm) 
    
    
    if (plot):
        im = plt.imshow(img, aspect="equal", origin="lower", cmap=matplotlib.cm.gray_r, interpolation="none", vmin=zmin, vmax=zmax)

        
        if (len(star_pix[:,0][mask]) >0):        
            plt.scatter(star_pix[:,0][mask], star_pix[:,1][mask], marker="o", s=np.minimum(150, 10000*(10./mag[mask][mask2])**9), edgecolor="red", facecolor="none")
        
        if (len(star_pix[:,0][mask][mask2]) >0):        
            plt.scatter(star_pix[:,0][mask][mask2], star_pix[:,1][mask][mask2], marker="o", s=20, edgecolor="yellow", facecolor="none")

        if (len(star_pix[:,0][mask][mask2][mask3]) >0):        
            plt.scatter(star_pix[:,0][mask][mask2][mask3], star_pix[:,1][mask][mask2][mask3], marker="o", s=200, edgecolor="green", facecolor="none")


        selected = star_pix[:,:][mask][mask2][mask3][mask_valid_fwhm]
        if (len(selected) >0):        
            plt.scatter(selected[:,0], selected[:,1], marker="o", \
                s=400, edgecolor="blue", facecolor="none")
            for i in np.arange(len(selected)):
                plt.text(selected[i,0]+10, selected[i,1]+10, i+1)
        
        plt.savefig(imfile.replace('.fits', '.seqstars.png'))
        print "Saved stars to ",imfile.replace('.fits', '.seqstars.png')
        plt.clf()
        
    return True
        
def add_to_zp_cal(ref_stars, image, logname):
    '''
    Records the parameters and instrumental and standard star magnitudes needed for zeropoint calibration.
    
    minst = M + c0 + c1(AIRMASS) + c2(color) + c3(UT)    
    '''
    
    coldic = {'u':'g', 'g':'r', 'r':'i', 'i':'z', 'z':'i'}

    r = np.genfromtxt(ref_stars, delimiter=" ", dtype=None, names=True)
    my = np.genfromtxt(image+".app.mag", comments="#", dtype=[("id","<f4"),  ("X","<f4"), ("Y","<f4"),("Xshift","<f4"), ("Yshift","<f4"),("fwhm","<f4"), ("ph_mag","<f4"), ("stdev","<f4"), ("fit_mag","<f4"), ("fiterr","<f4")])

    if (my.size < 2):
        my = np.array([my])
    if (r.size < 2):
        r = np.array([r])
    coldic = {'u':'g', 'g':'r', 'r':'i', 'i':'z', 'z':'i'}
    
    band = fitsutils.get_par(image, 'filter')
    mask_valid1 = np.array(my["fwhm"]<9000) * np.array(my["ph_mag"]<9000) * np.array(~ np.isnan(r[band])) * np.array(~ np.isnan(my["fit_mag"]))
    
    r = r[mask_valid1]
    my = my[mask_valid1]
    N = len(r)
    
    my["fiterr"][np.isnan( my["fiterr"])] = 100

    col_band = coldic[band]
    exptime = fitsutils.get_par(image, "exptime")
    date = fitsutils.get_par(image, "JD")
    airmass = fitsutils.get_par(image, "AIRMASS")
    
    if (not os.path.isfile(logname)):
            with open( logname, "a") as f:
                f.write("#filename,filter,std,stderr,inst,insterr,jd,airmass,color,exptime\n")
    with open( logname, "a") as f:
        for i in range(N):
            f.write("%s,%s,%.3f,%.3f,%3f,%.3f,%.3f,%.3f,%.3f,%.3f\n"%
            (image, band, r[i][band], r[i]["d"+band], my[i]['fit_mag'], my[i]['fiterr'], date, airmass, r[i][band]-r[i][col_band],exptime))
       
def lsq_zeropoint(logfile, plot=True):
    '''
    Uses least squares approach to compute the optimum coefficients for ZP, colour term, airmass and time.
    
    '''
    a = np.genfromtxt(logfile, dtype=None, names=True, delimiter=",")
    
    a['jd'] = a['jd'] - np.min(a['jd'])
    cols = {'u':'purple', 'g':'green', 'r':'red', 'i':'orange'}


    for b in set(a['filter']):
        ab = a[a['filter']==b]
        M = np.zeros((len(ab), 4))
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
        
        if (plot):
            f, ((ax1, ax2), (ax3, ax4)) = plt.subplots(2, 2)
            ax1.plot(ab['color'], ab['std'] - ab['inst'], "o", color=cols[b])
            ax1.plot(ab['color'], coef[0] + ab['color']*coef[1], color=cols[b])
            ax1.set_xlabel("color")

            ax2.plot(ab['airmass'], ab['std'] - ab['inst'], "o", color=cols[b])
            ax2.plot(ab['airmass'], coef[0] + (1.3+ab['airmass'])*coef[2], color=cols[b])
            ax2.set_xlabel("airmass")

            ax3.plot(ab['jd'], ab['std'] - ab['inst'], "o", color=cols[b])
            ax3.plot(ab['jd'], coef[0] + ab['jd']*coef[3], color=cols[b])
            ax3.set_xlabel("jd")

            
    if (plot):
        plt.show()
def find_zeropoint_noid(ref_stars, image, plot=False, psf=False):
    '''
    Finds the zeropoint by comparig the magnitude of the stars measured in the image,
    vs. the magnitude of the stars from the reference catalogue.
    band: band for which we wat to measure the zeropoint.
    col_band: band for the colour term we want to use.
    psf: points whether we want to use aperture of psf photometry on sequence stars.
    
    '''
    
    print 'Finding the optimum ZP fit...'
    
    r = np.genfromtxt(ref_stars, delimiter=" ", dtype=None, names=True)
    my = np.genfromtxt(image+".app.mag", comments="#", dtype=[("id","<f4"),  ("X","<f4"), ("Y","<f4"),("Xshift","<f4"), ("Yshift","<f4"),("fwhm","<f4"), ("ph_mag","<f4"), ("stdev","<f4"), ("fit_mag","<f4"), ("fiterr","<f4")])

    if (my.size < 2):
        my = np.array([my])
    if (r.size < 2):
        r = np.array([r])
    coldic = {'u':'g', 'g':'r', 'r':'i', 'i':'z', 'z':'i'}
    band = fitsutils.get_par(image, 'filter')
    col_band = coldic[band]
    
    mask_valid1 = np.array(my["fwhm"]<9000) * np.array(my["ph_mag"]<9000) * np.array(~ np.isnan(r[band])) * np.array(~ np.isnan(my["fit_mag"]))
    
    N = len(r)
    r = r[mask_valid1]
    my = my[mask_valid1]
    
    #mask_valid1 = mask_valid1 * np.abs(r[band]-r[col_band])<1
    #r = r[mask_valid1]
    #my = my[mask_valid1]


    my["fiterr"][np.isnan( my["fiterr"])] = 100
    
    ids = np.arange(N)+1
    ids = ids[mask_valid1]
    #print ids
    #print r[band]-r[col_band], r[band] - my["fit_mag"], my["fit_mag"]

    coefs, residuals, rank, singular_values, rcond = np.polyfit(r[band]-r[col_band], r[band] - my["fit_mag"], w=1./np.maximum(0.1, np.sqrt(my["fiterr"]**2 + r['d'+band]**2)), deg=1, full=True)
    p = np.poly1d(coefs)
    
    if (plot):
        print "Plotting..."
        plt.errorbar(r[band]-r[col_band], r[band] - my["fit_mag"] , yerr=np.sqrt(my["fiterr"]**2 + r['d'+band]**2), marker="o", ls="None")
        for i, myid in enumerate(ids):
            plt.text(r[band][i]-r[col_band][i] + 0.01, r[band][i] - my["fit_mag"][i]+ 0.01, str(myid))
        x = np.linspace(np.min(r[band]-r[col_band]), np.max(r[band]-r[col_band]), 100)
        plt.plot(x, p(x))
        plt.title("Best fit ZP:"+str( p[0])+ " colour term " + str(p[1]))
        plt.xlabel("{:} - {:}".format(band, col_band))
        plt.ylabel("ZP")
        plt.savefig(image.replace(".fits", ".zp.png"))
        plt.clf()
    
    print band, "-", col_band, p[0], p[1]
    
    pred = p(r[band]-r[col_band])
    return p[0], p[1], np.sqrt(np.sum( (pred - (r[band]- my["fit_mag"]))**2))/len(pred)

def calibrate_zeropoint(image, plot=True, debug=False, refstars=None):
    
    filt = fitsutils.get_par(image, 'filter')
    exptime = fitsutils.get_par(image, "exptime")
    if fitsutils.has_par(image, "JD"):
        date = fitsutils.get_par(image, "JD")
    elif fitsutils.has_par(image, "MJD"):
        date = fitsutils.get_par(image, "MJD")

    objname = fitsutils.get_par(image, "OBJECT")
    airmass = fitsutils.get_par(image, "AIRMASS")
    
    print "Starting calibration of ZP for image", image,"for object", objname

    if (exptime < 10):
        print "ERROR. Exposure time too short for this image to see anything..."
        return 

    extracted = extract_star_sequence(image, filt, plot=plot, survey='sdss', debug=debug, refstars=refstars)
    if (not extracted):
        print "Field not in SDSS or error when retrieving the catalogue... Skipping."
        return 

    #If extraction worked, we can get the FWHM        
    fwhm = fitsutils.get_par(image, "fwhm")
    fwhm_as = fwhm * 0.394

    app_phot.get_app_phot("/tmp/sdss_cat_det.txt", image, wcsin='logic')
    z, c, err = find_zeropoint_noid("/tmp/sdss_cat_det.txt", image, plot=plot)
    
    #Log the current zeropoint for this image
    logname = os.path.join(os.path.dirname(image), "zeropoint.log")
    
    #Add the data to a later stage zeropoint calibrtion with all-sky data.
    zplogname = os.path.join(os.path.dirname(image), "allstars_zp.log")
    
    add_to_zp_cal("/tmp/sdss_cat_det.txt", image, zplogname)


    if (not os.path.isfile(logname)):
            with open( logname, "a") as f:
                f.write("#filename,exptime,filter,date,airmass,fwhm_pix,fwhm_as,zeropoint,color,err\n")
    with open( logname, "a") as f:
        f.write("%s,%.1f,%s,%3f,%.3f,%.3f,%.3f,%.3f,%.3f,%.3f\n"%(image,exptime,filt,date,airmass,fwhm,fwhm_as,z,c,err))
        
def plot_zp(zpfile):
    import datetime
    
    def mkdate(text):
        return datetime.datetime.strptime(text, '%Y-%m-%dT%H:%M:%S') 
    
    a = np.genfromtxt(zpfile, names=True, dtype=None, delimiter=',')
    
    cols = {'u':'purple', 'g':'green', 'r':'red', 'i':'orange'}
    
    for fi in ['u', 'g', 'r', 'i']:
        for i in range(len(a[a['filter']==fi])):
            plt.errorbar( (a[a['filter']==fi]['date'][i]-np.min(a[a['filter']==fi]['date'], axis=0))*24., \
            a[a['filter']==fi]['zeropoint'][i], yerr=a[a['filter']==fi]['err'][i], marker='o', mfc=cols[fi], mec='k', ecolor=cols[fi], ls='none', ms=a[a['filter']==fi]['fwhm_pix'][i])
        print "Median zeropoint for filter %s: %.2f mag"%(fi, np.median(a[a['filter']==fi]['zeropoint']))

    plt.gca().invert_yaxis()
    plt.xlabel("Obs Date (JD - min(JD)) [h]")
    plt.ylabel("ZP [mag]")
    plt.show()
    
def plot_zp_airmass(zpfile):
    import datetime
    
    def mkdate(text):
        return datetime.datetime.strptime(text, '%Y-%m-%dT%H:%M:%S') 
    
    a = np.genfromtxt(zpfile, names=True, dtype=None, delimiter=',')
    
    cols = {'u':'purple', 'g':'green', 'r':'red', 'i':'orange'}
    
    for fi in ['u', 'g', 'r', 'i']:
        for i in range(len(a[a['filter']==fi])):
            plt.errorbar( a[a['filter']==fi]['airmass'][i], a[a['filter']==fi]['zeropoint'][i], yerr=a[a['filter']==fi]['err'][i], marker='o', mfc=cols[fi], mec='k', ecolor=cols[fi], ls='none', ms=a[a['filter']==fi]['fwhm_pix'][i])
        print "Median zeropoint for filter %s: %.2f mag"%(fi, np.median(a[a['filter']==fi]['zeropoint']))

    plt.gca().invert_yaxis()
    plt.xlabel("Airmass")
    plt.ylabel("ZP [mag]")
    plt.show()