# -*- coding: utf-8 -*-
"""
Created on Tue Jul 14 15:01:50 2015

@author: nadiablago
"""
import pyfits as pf
import zscale
import pywcs
import numpy as np
import matplotlib
from matplotlib import pylab as plt
from scipy import ndimage
import aplpy
import coordinates_conversor
import fitsmanip
import fitsutils
import time_utils

    
def finder(myfile,searchrad=0.2/60.):
    
    ra, dec = coordinates_conversor.hour2deg(fitsutils.get_par(myfile, "RA"), fitsutils.get_par(myfile, "DEC"))


    hdulist = pf.open(myfile)[0]
    img = hdulist.data * 1.            
    
    wcs = pywcs.WCS(hdulist.header)

    target_pix = wcs.wcs_sky2pix([(np.array([ra,dec], np.float_))], 1)[0]
    corner_pix = wcs.wcs_sky2pix([(np.array([ra,dec+searchrad], np.float_))], 1)[0]
    dx = int(np.abs(np.ceil(corner_pix[1] - target_pix[1])))
    
    imgslice = img[target_pix[0]-2*dx:target_pix[0]+2*dx, target_pix[1]-2*dx:target_pix[1]+2*dx]
    #zmin, zmax = zscale.zscale()
    zmin = np.percentile(imgslice.flatten(), 5)
    zmax = np.percentile(imgslice.flatten(), 99.9)
    
    gc = aplpy.FITSFigure(myfile, figsize=(10,9), north=True)
    gc.show_grayscale(vmin=zmin, vmax=zmax)
    gc.show_scalebar(0.1/60.)
    gc.scalebar.set_label('10 arcsec')
    gc.scalebar.set_color('white')
    gc.recenter(ra, dec, searchrad)
    gc.show_markers(ra,dec+searchrad/20.,edgecolor='red',facecolor='none',marker="|",s=250, lw=4)
    gc.show_markers(ra-(searchrad/20.)/np.cos(np.deg2rad(dec)),dec,edgecolor='red',facecolor='none',marker="_",s=250, lw=4)

    ras = np.array([ra+searchrad*0.7/ np.cos(np.deg2rad(dec)), ra+searchrad*0.7/ np.cos(np.deg2rad(dec))])
    decs = np.array([dec-searchrad*0.9, dec-searchrad*0.9])
    dxs = np.array([0, searchrad/5 / np.cos(np.deg2rad(dec))])
    dys = np.array([searchrad/5, 0])

    gc.show_arrows(ras, decs, dxs, dys, edgecolor="k", facecolor="k")
    gc.add_label(ras[0]+dxs[0]*1.1, decs[0]+dys[0]*1.1, 'N', relative=False, color="k", horizontalalignment="center")
    gc.add_label(ras[1]+dxs[1]*1.1, decs[1]+dys[1]*1.1, 'E', relative=False, color="k", horizontalalignment="center")


    name = fitsutils.get_par(myfile, "OBJECT")
    gc.add_label(0.05, 0.95, 'Object: %s'%(name), relative=True, color="white", horizontalalignment="left")                   
    gc.add_label(0.05, 0.9, 'Coordinates: RA=%s DEC=%s'%(coordinates_conversor.deg2hour(ra, dec)), relative=True, color="white", horizontalalignment="left")
    gc.add_label(0.05, 0.84, 'Filter: SDSS r', relative=True, color="white", horizontalalignment="left")
    
    gc.save(myfile.replace(".fits", '_%s_finder.jpg'%(name)))
    
