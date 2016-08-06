import glob
import os
import time
import numpy as np


def report():
    """Generate DRP report using output sp_<object>.npy files"""

    flist = [f for f in glob.glob("sp_*.npy")
             if "_A_" not in f and "_B_" not in f]
    print "\nReport generated on %s" % time.strftime("%c")
    print "\nSEDM DRP run in %s\nFound %d sp_*.npy files\n" % \
          (os.getcwd(), len(flist))
    totexpt = 0.
    print "Object                     Obs  Exptime Qual Skysb"
    for f in flist:
        if '_A_' in f or '_B_' in f:
            continue
        sp = np.load(f)[0]
        if 'quality' in sp:
            qual = sp['quality']
        else:
            qual = 0
        if 'sky_subtraction' in sp:
            skysub = sp['sky_subtraction']
        else:
            skysub = 1
        if '_obs' in f:
            obs = f.split('_')[2].split('.')[0]
        else:
            obs = "A/B"

        if 'exptime' in sp:
            expt = sp['exptime']
            if "A/B" in obs:
                expt *= 2.
        else:
            expt = 0.
        totexpt += expt

        objname = f.split('_')[1].split('.')[0]

        print "%-25s %4s   %6.1f %4d %5s" % (objname, obs, expt, qual,
                                           ("on" if skysub else "off"))
    print "\nTotal science exposure time = %.1f s" % totexpt

if __name__ == '__main__':
    report()