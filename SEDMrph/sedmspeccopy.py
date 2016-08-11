# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:23:00 2016

@author: nblago
"""
import datetime
import glob, os
import argparse
import subprocess
import re



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''

        Runs astrometry.net on the image specified as a parameter and returns 
        the offset needed to be applied in order to center the object coordinates 
        in the reference pixel.
            
        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('-d', '--photdir', type=str, dest="photdir", help='Fits directory file with tonight images.', default=None)

    args = parser.parse_args()
    
    photdir = args.photdir
    print "Parameter directory where stats are run :",photdir
    
    if (photdir is None):
        timestamp=datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-","")
        photdir = os.path.join("/scr2/sedmdrp/redux/", timestamp)

    os.chdir(photdir)
    
    sedmfiles = glob.glob("PTF*.txt")
    print "Copying", sedmfiles
    
    for f in sedmfiles:
        
        qual = -1
        #retrieve the quality of the spectra.
        with open(f, "r") as sf:
            a = sf.readlines()
            
            qual = [ai for ai in a if "QUALITY" in ai]
            
            if len(qual) > 0:
                match = re.search(r'\(?([0-9]+)\)?', qual[0])
                qual = int(match.group(1))

        #Only write the spectra that have good qualities.
        if (qual < 3):
            newname = os.path.basename(f).split("_")[0].replace("PTF", "")
            newname += "_%s_P60_v1.ascii"%timestamp
            cmd = "rcp %s sedm@yupana.caltech.edu:/scr/apache/htdocs/marshals/transient/ptf/spectra/sedm_to_upload/%s"%(f, newname)
            print cmd
            try:
                subprocess.call(cmd, shell=True)
            except:
                print "Error copying the file"
                pass
