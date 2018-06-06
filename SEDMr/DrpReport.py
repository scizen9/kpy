import glob
import os
import time
import numpy as np


def report():
    """Generate DRP report using output sp_<object>.npy files"""

    flist = [f for f in glob.glob("sp_*.npy")
             if "_A_" not in f and "_B_" not in f]
    flist.sort()
    print("\nReport generated on %s" % time.strftime("%c"))
    print("\nSEDM DRP run in %s\nFound %d sp_*.npy files\n" %
          (os.getcwd(), len(flist)))
    print("\nSee http://pharos.caltech.edu/data_access/ifu?date=%s" %
          os.getcwd().split('/')[-1])
    totexpt = 0.
    lostexp = 0.
    print("Object                     Obs Method  Exptime Qual Skysb Airmass "
          "   Reducer   Type      z         Rlap")
    for f in flist:
        if '_A_' in f or '_B_' in f:
            continue
        # trim the .npy off the end
        objname = '.'.join(f.split('.')[0:-1])
        # check the ascii spectrum file for a type
        sfile = objname[3:] + "_SEDM.txt"
        ctype = ""
        zmch = ""
        rlap = ""
        if os.path.exists(sfile):
            with open(sfile, "r") as sfl:
                lines = sfl.readlines()
                # check for previous classification
                clas = [li for li in lines if "TYPE" in li]
                ctype = ""
                if len(clas) > 0:
                    for cl in clas:
                        ctype += (" %s" % cl.split()[-1])
                # get redshift
                zmch = [li for li in lines if "REDSHIFT" in li]
                if len(zmch) > 0:
                    zmch = ("%.4f" % float(zmch[0].split()[-1]))
                else:
                    zmch = ""
                # get rlap
                rlap = [li for li in lines if "RLAP" in li]
                if len(rlap) > 0:
                    rlap = ("%.2f" % float(rlap[0].split()[-1]))
                else:
                    rlap = ""
                sfl.close()

        # load the spectrum file
        sp = np.load(f)[0]
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

        print("%-25s %4s %6s   %6.1f %4d %5s  %5.3f   %9s  %-9s  %6s  %6s" %
              (objname, obs, meth, expt, qual, ("on" if skysub else "off"),
               air, reducer, ctype, zmch, rlap))
    print("\nTotal quality (1-3) science exposure time = %.1f s" % totexpt)
    if lostexp > 0:
        print("Total exposure time lost to bad targets = %.1f s\n" % lostexp)


if __name__ == '__main__':
    report()
