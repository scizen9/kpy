import glob
import os
import time
import numpy as np


def report():
    """Generate DRP report using output sp_<object>.npy files"""

    flist = [f for f in glob.glob("sp_*.npy")
             if "_A_" not in f and "_B_" not in f]
    flist.sort(key=os.path.getmtime)
    print("\nReport generated on %s" % time.strftime("%c"))
    print("\nSEDM DRP run in %s\nFound %d sp_*.npy files\n" %
          (os.getcwd(), len(flist)))
    totexpt = 0.
    lostexp = 0.
    print("Object                     Obs Method  Exptime Qual Skysb Airmass "
          "Reducer")
    for f in flist:
        if '_A_' in f or '_B_' in f:
            continue
        sp = np.load(f)[0]
        # trim the .npy off the end
        objname = '.'.join(f.split('.')[0:-1])
        if 'quality' in sp:
            qual = sp['quality']
        else:
            qual = 0
        if 'reducer' in sp:
            reducer = sp['reducer']
        else:
            reducer = '-'
        if 'sky_subtraction' in sp:
            skysub = sp['sky_subtraction']
        else:
            skysub = 1
        if '_obs' in f:
            if len(objname.split('_')) > 2:
                obs = objname.split('_')[-1]
            elif len(objname.split('_')) == 2:
                obs = objname.split('_')[-1]
            else:
                obs = "obs1"
        else:
            obs = "obs1"
        obs = obs.split('s')[-1]
        if 'object_spaxel_ids_A' in sp:
            meth = "A / B"
        else:
            meth = "Single"

        if 'exptime' in sp:
            expt = sp['exptime']
            if "A / B" in meth:
                expt *= 2.
        else:
            expt = 0.
        # get airmass
        meta = sp['meta']
        if 'airmass1' in meta:
            air = meta['airmass1']
            if 'airmass2' in meta:
                air = (air + meta['airmass2']) / 2.
        elif 'airmass' in meta:
            air = meta['airmass']
        else:
            air = 0.

        # Don't count bad objects
        if qual < 3:
            totexpt += expt
        else:
            lostexp += expt

        if '_obs' in objname:
            objname = "_".join(objname.split('_')[1:-1])
        else:
            objname = "_".join(objname.split('_')[1:])

        print("%-25s %4s %6s   %6.1f %4d %5s  %5.3f   %s" % (objname, obs, meth,
                                                             expt, qual,
                                                             ("on" if skysub
                                                              else "off"), air,
                                                             reducer))
    print("\nTotal quality (1-3) science exposure time = %.1f s" % totexpt)
    if lostexp > 0:
        print("Total exposure time lost to bad targets = %.1f s\n" % lostexp)

if __name__ == '__main__':
    report()
