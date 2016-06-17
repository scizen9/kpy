# -*- coding: utf-8 -*-
"""
Created on Wed Jun 15 20:31:17 2016

@author: nadiablago
"""

import sextractor
import glob, os
import numpy as np
from astropy.coordinates import SkyCoord
from matplotlib import pylab as plt
import argparse
import pyfits as pf


plotdir="/Users/nadiablago/Documents/Projects/SEDM/flex"

def get_par(myfits, par):
    '''
    Returns the header parameter from the fits.
    '''
    hdu = pf.open(myfits)
    header = hdu[0].header
    return header[str.upper(par)]
    
def run_flexure_test(sexfiles, plotdir=plotdir):
    
    posfiles = []
    
    for sf in sexfiles:
        c = np.genfromtxt(sf)
        
        # We want bright sources.
        c = c[c[:,2]<-10]
        
        # We want round-shaped sources.
        c = c[c[:,7]<0.6]
        
        #Save the positions in a separate file with only X Y.
        np.savetxt(sf.replace(".sex", ".pos"), c[:,0:2])
        posfiles.append(sf.replace(".sex", ".pos"))
        
    c0 = np.genfromtxt(posfiles[0])
    for f in posfiles[1:]:
        print f
        c1 = np.genfromtxt(f)
    
        c = SkyCoord(x=c0[:,0], y=c0[:,1], z=np.zeros(len(c0)), unit='m', representation='cartesian')
        catalog = SkyCoord(x=c1[:,0], y=c1[:,1], z=np.zeros(len(c1)), unit='m', representation='cartesian')#SkyCoord(x=c1[:0]*u.pixel, y=c1[:,1]*u.pixel, z=np.zeros(len(c0))*u.pixel, representation='cartesian')  
    
        #Now idx are indices into catalog that are the closest objects to each of the 
        #coordinates in c, d2d are the on-sky distances between them.
        idx, d2d, d3d = c.match_to_catalog_sky(catalog)  
        
        matches = catalog[idx]
        matches = matches[d3d.value<10]
        d3d = d3d[d3d.value<10]

        plt.hist(d3d, bins=50)
        plt.xlabel("Deviation [pixels]")
        plt.title("Median deviation: %.3f"%np.median(d3d.value))
        plt.savefig(os.path.join(plotdir, "%s_vs_%s.png"%(os.path.basename(posfiles[0]), os.path.basename(f))))
        plt.clf()

        plt.scatter(matches.x.value, matches.y.value, c=np.minimum(np.mean(d3d.value)+2*np.std(d3d.value), d3d.value))
        c = plt.colorbar()
        c.set_label("Deviation [pixels]")
        plt.savefig(os.path.join(plotdir, "%s_vs_%s_XY.png"%(os.path.basename(posfiles[0]), os.path.basename(f))))
        plt.clf()

     
if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description=\
        '''

        Runs astrometry.net on the image specified as a parameter and returns 
        the offset needed to be applied in order to center the object coordinates 
        in the reference pixel.
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('raw', type=str, help='Directory containing the raw fits for the night.')

    args = parser.parse_args()
    
    raw = args.raw
    
    if (raw is None):
        print "Please, add the directory containing raw data as a parameter."
    else:
        files = glob.glob(os.path.join(raw, "ifu*fits"))
        files_hg = [f for f in files if "Calib:  Hg" in get_par(f, "OBJECT")]
        
        if (len(files_hg)>1):
            sexfiles = sextractor.run_sex(files_hg, mask=False)
            
            plotdir = os.path.join(raw, "stats")
            if (not os.path.isdir(plotdir)):
                os.makedirs(plotdir)
                
            run_flexure_test(sexfiles, plotdir=plotdir)
