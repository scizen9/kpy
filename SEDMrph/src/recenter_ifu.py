# -*- coding: utf-8 -*-
"""
Created on Thu Jan  7 11:30:32 2016

@author: nadiablago
"""
import fitsutils
import subprocess, os, sys
import pyfits as pf
import pywcs
import coordinates_conversor as cc
import numpy as np
import argparse

def solve_astrometry(img, radius=6.5, with_pix=True, first_call=True):
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
        cmd = cmd + " --scale-units arcsecperpix  --scale-low 0.375 --scale-high 0.425"
    print cmd
    
    cmd = cmd + " > /dev/null"

    subprocess.call(cmd, shell=True)
    
    #Cleaning after astrometry.net
    if (os.path.isfile(img.replace(".fits", ".axy"))):
        os.remove(img.replace(".fits", ".axy"))
    if (os.path.isfile(img.replace(".fits", "-indx.xyls"))):
        os.remove(img.replace(".fits", "-indx.xyls"))
    if (os.path.isfile("none")):
        os.remove("none")
        
        
    if(not os.path.isfile(astro) and first_call):
        print "Astrometry failed on file %s! Trying with a larger radius..."%astro
        solve_astrometry(img, radius=8, with_pix=True, first_call=False)
        
    return astro
    
def get_offset_ref_pixel(f):
    
    if(not os.path.isfile(f)):
        print "File %s does not exist! Returning Zero offsets..."%f
        return 0,0
    else:
        image = pf.open(f)
        wcs = pywcs.WCS(image[0].header)
        rra, rdec = cc.hour2deg(image[0].header['OBRA'],image[0].header['OBDEC'] )
        pra, pdec = wcs.wcs_pix2sky(np.array([[1293., 1280.]] , np.float_), 1)[0]
        dra, ddec = cc.get_offset(pra, pdec, rra, rdec)
        
        return dra, ddec

def main(infile):
    
    offset_file = "/tmp/%s_dra_ddec.txt"%(os.path.basename(infile).replace(".fits", ""))

    if (os.path.isfile(offset_file)):
        print "Offset file %s already exists for the image."%offset_file
        return
        
    newfile = solve_astrometry(infile)
    dra, ddec = get_offset_ref_pixel(newfile)
    
    
    np.savetxt(offset_file, np.array([("A", dra-10, ddec), ("B", dra+5, ddec-2)]), fmt="%s")
    
    print "Offsets computed: \n A %.4f %.4f \n B %.4f %.4f"%(dra-10, ddec, dra+5, ddec-2)
    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

        Runs astrometry.net on the image specified as a parameter and returns 
        the offset needed to be applied in order to center the object coordinates 
        in the reference pixel.
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('image', type=str, help='Fits file with acquisition image.')

    args = parser.parse_args()
    
    infile = args.image
    
    main(infile)
