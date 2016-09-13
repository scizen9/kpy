# -*- coding: utf-8 -*-
"""
Created on Thu Mar  3 19:23:29 2016

@author: nadiablago
"""

import matplotlib
matplotlib.use("Agg")
import rcred
import subprocess
import glob, os, time
import argparse
import fitsutils
import datetime
import zeropoint
import logging

#Log into a file
FORMAT = '%(asctime)-15s %(levelname)s [%(name)s] %(message)s'
#root_dir = "/tmp/logs/"
root_dir = "/scr2/sedm/logs/"
now = datetime.datetime.utcnow()
timestamp=datetime.datetime.isoformat(now)
timestamp=timestamp.split("T")[0]
logging.basicConfig(format=FORMAT, filename=os.path.join(root_dir, "rcred_{0}.log".format(timestamp)), level=logging.INFO)
logger = logging.getLogger('realtimered')


def reduce_all_dir(photdir, overwrite=False):
    
    #Reduce the data that is already in the directory.
    cmd = "python %s/rcred.py -d %s"%(os.environ["SEDMPH"], photdir)    
    if (overwrite):
        cmd = cmd + " -o"
    subprocess.call(cmd, shell=True)
    logger.info("Reduce all dir: %s"%cmd)

    # Copy the content of the reduced directory into a new directory with the date of the observations.
    dayname = os.path.basename(photdir)
    reducedname = os.path.join(photdir, "reduced")

    #Reduce the data that is already in the directory.
    cmd = "python %s/zeropoint.py  %s"%(os.environ["SEDMPH"], reducedname)    
    subprocess.call(cmd, shell=True)
    logger.info("zeropoint for all dir: %s"%cmd)
    
    if (os.path.isdir(reducedname)):
    	cmd = "rcp -r %s grbuser@transient.caltech.edu:/scr3/mansi/ptf/p60phot/fremling_pipeline/sedm/reduced/%s"%(reducedname, dayname)
    	subprocess.call(cmd, shell=True)
    else:
	os.makedirs(reducedname)
    
def reduce_on_the_fly(photdir):
    '''
    Waits for new images to appear in the directory to trigger their incremental reduction as well.
    '''
    #Get the current the number of files
    nfiles = glob.glob(os.path.join(photdir, "rc*[0-9].fits"))
    
    dayname = os.path.basename(photdir)
    
    time_ini = datetime.datetime.now()
    time_curr = datetime.datetime.now()
    
    #Run this loop for 12h since the start.
    while (time_curr-time_ini).total_seconds() < 12.*3600:
        nfilesnew = glob.glob(os.path.join(photdir, "rc*[0-9].fits"))
        if len(nfilesnew) == len(nfiles):
            time.sleep(10)
        else:
            new = [f for f in nfilesnew if f not in nfiles]
            for n in new:
                if (not fitsutils.has_par(n, "IMGTYPE")):
                    print "Image",n,"Does not have an IMGTYPE"
                    time.sleep(0.5)
                    if (not fitsutils.has_par(n, "IMGTYPE")):
                        print "Image",n,"STILL Does not have an IMGTYPE"
                        continue
                if (fitsutils.get_par(n, "IMGTYPE")=="SCIENCE"):
                    reduced = rcred.reduce_image(n)
                    try:
                        zeropoint.calibrate_zeropoint(reduced)
                    except:
                        logger.error("Could not calibrate zeropoint for image %s"%reduced)
                    #Copy them to transient
                    for r in reduced:
                        cmd = "rcp %s grbuser@transient.caltech.edu:/scr3/mansi/ptf/p60phot/fremling_pipeline/sedm/reduced/%s/."%(r, dayname)
                        subprocess.call(cmd, shell=True)
                        logger.info(cmd)
                        logger.info( "Successfully copied the image: %s"% cmd)
        time_curr = datetime.datetime.now()
        nfiles = nfilesnew  
        
    
         
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

        Runs astrometry.net on the image specified as a parameter and returns 
        the offset needed to be applied in order to center the object coordinates 
        in the reference pixel.
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-d', '--photdir', type=str, dest="photdir", help='Fits directory file with tonight images.', default=None)
    parser.add_argument('-f', '--fullred', action='store_true', dest="fullred", default=False, help='Whether we should do a full reduction.')
    parser.add_argument('-o', '--overwrite', action="store_true", help='re-reduce and overwrite the reduced images?', default=False)

    args = parser.parse_args()
    
    photdir = args.photdir
    fullred = args.fullred
    overwrite = args.overwrite
    
    if (photdir is None):
        timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-","")
        photdir = os.path.join("/scr2/sedm/phot/", timestamp)
    if (fullred):
        reduce_all_dir(os.path.abspath(photdir), overwrite=overwrite)
    reduce_on_the_fly(os.path.abspath(photdir))
    
    #After 12h, invoke the zeropoint calibration.
    zeropoint.main(os.path.join(os.path.abspath(photdir), "reduced"))
