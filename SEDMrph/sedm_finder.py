# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 15:01:50 2015

@author: nadiablago
"""
import matplotlib
matplotlib.use("Agg")
import pyfits as pf
import zscale
import pywcs
import numpy as np
import aplpy
import coordinates_conversor
import fitsutils
import datetime
import os
import glob
import argparse
import subprocess
from matplotlib import pylab as plt
 
def finder(myfile,searchrad=0.2/60.):
    
    ra, dec = coordinates_conversor.hour2deg(fitsutils.get_par(myfile, "OBJRA"), fitsutils.get_par(myfile, "OBJDEC"))


    hdulist = pf.open(myfile)[0]
    img = hdulist.data * 1.            
    img = img.T

    wcs = pywcs.WCS(hdulist.header)

    target_pix = wcs.wcs_sky2pix([(np.array([ra,dec], np.float_))], 1)[0]
    corner_pix = wcs.wcs_sky2pix([(np.array([ra,dec+searchrad], np.float_))], 1)[0]
    dx = int(np.abs(np.ceil(corner_pix[1] - target_pix[1])))
    
    imgslice = img[int(target_pix[0])-2*dx:int(target_pix[0])+2*dx, int(target_pix[1])-2*dx:int(target_pix[1])+2*dx]
    #zmin, zmax = zscale.zscale()
    zmin = np.percentile(imgslice.flatten(), 5)
    zmax = np.percentile(imgslice.flatten(), 98)
   
    print "Min: %.1f, max: %.1f"%(zmin, zmax) 
    gc = aplpy.FITSFigure(myfile, figsize=(10,9), north=True)
    gc.show_grayscale(vmin=zmin, vmax=zmax, smooth=1, kernel="gauss")
    gc.show_scalebar(0.1/60.)
    gc.scalebar.set_label('10 arcsec')
    gc.scalebar.set_color('white')
    gc.recenter(ra, dec, searchrad)
    #gc.show_markers(ra,dec+searchrad/20.,edgecolor='red',facecolor='none',marker="|",s=250, lw=10)
    #gc.show_markers(ra-(searchrad/20.)/np.cos(np.deg2rad(dec)),dec,edgecolor='red',facecolor='none',marker="_",s=250, lw=10)

    ras = np.array([ra , ra])
    decs = np.array([dec, dec])
    dxs = np.array([0, searchrad/10 / np.cos(np.deg2rad(dec))])
    dys = np.array([searchrad/10, 0])

    gc.show_arrows(ras, decs, dxs, dys, edgecolor="red", facecolor="red", head_width=0)

    ras = np.array([ra+searchrad*0.7/ np.cos(np.deg2rad(dec)), ra+searchrad*0.7/ np.cos(np.deg2rad(dec))])
    decs = np.array([dec-searchrad*0.9, dec-searchrad*0.9])
    dxs = np.array([0, searchrad/5 / np.cos(np.deg2rad(dec))])
    dys = np.array([searchrad/5, 0])

    gc.show_arrows(ras, decs, dxs, dys, edgecolor="k", facecolor="k")
    gc.add_label(ras[0]+dxs[0]*1.1, decs[0]+dys[0]*1.1, 'N', relative=False, color="k", horizontalalignment="center")
    gc.add_label(ras[1]+dxs[1]*1.1, decs[1]+dys[1]*1.1, 'E', relative=False, color="k", horizontalalignment="center")


    name = fitsutils.get_par(myfile, "NAME")
    filter = fitsutils.get_par(myfile, "FILTER")
    gc.add_label(0.05, 0.95, 'Object: %s'%(name), relative=True, color="white", horizontalalignment="left")                   
    gc.add_label(0.05, 0.9, 'Coordinates: RA=%s DEC=%s'%(coordinates_conversor.deg2hour(ra, dec)), relative=True, color="white", horizontalalignment="left")
    gc.add_label(0.05, 0.84, 'Filter: SDSS %s'%filter, relative=True, color="white", horizontalalignment="left")
    
    findername = 'finder_%s_%s.jpg'%(name, filter)
    gc.save(findername)
    
    return findername
    

if __name__=="__main__":  
    parser = argparse.ArgumentParser(description=\
    '''

    Creates a finder chart for every acquisition image in the folder specified as a parameter.
    As a final step, it copies the acquisition image to the "agn" machine to visualize it.
        
    ''', formatter_class=argparse.RawTextHelpFormatter)
    
    parser.add_argument('-d', '--photdir', type=str, dest="photdir", help='Fits directory file with tonight images.', default=None)

    args = parser.parse_args()
    
    photdir = args.photdir
    
    
    if (photdir is None):
    	timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
    	timestamp = timestamp.split("T")[0].replace("-","")
        photdir = os.path.join("/scr2/sedm/phot/", timestamp)
    else:
	timestamp = os.path.basename(os.path.abspath(photdir))

    os.chdir(photdir)
    print "Changed to directory where the data is. %s"%photdir
    
    #We only generate onle finder with the first image.
    files = glob.glob("a_*fits")
    files.sort()
    for f in files:
	object = fitsutils.get_par(f, "OBJECT")
        if (fitsutils.get_par(f, "IMGTYPE") == "ACQUISITION" or fitsutils.get_par(f, "IMGTYPE") == "SCIENCE" ) and "STD" not in object and "BD" not in object and "SA" not in object:
	    findername = 'finder_%s_%s.jpg'%(fitsutils.get_par(f, "NAME"), fitsutils.get_par(f, "FILTER"))
	    if(not os.path.isfile(findername)):
            	findername = finder(f)
            if(os.path.isfile(findername)):
                cmd = "rcp %s sedm@agn.caltech.edu:/usr/apache/htdocs/astro/sedm/stats/%s/."%(findername, timestamp)
		print cmd
		subprocess.call(cmd, shell=True)

    cmd = "rcp /tmp/finders.php sedm@agn.caltech.edu:/usr/apache/htdocs/astro/sedm/stats/%s/."%(timestamp)
    print cmd
    subprocess.call(cmd, shell=True)

