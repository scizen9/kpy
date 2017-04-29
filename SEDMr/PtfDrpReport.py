import glob
import os
import time
import numpy as np
import subprocess


def report():
    """Generate DRP report using output sp_<object>.npy files"""

    flist = [f for f in glob.glob("sp_*.npy")
             if "_A_" not in f and "_B_" not in f]
    flist.sort(key=os.path.getmtime)
    
    out = open("report_ptf.txt", "w")
    out.write("\nReport generated on %s\n\n" % time.strftime("%c"))
    totexpt = 0.
    lostexp = 0.
    out.write("Object                     Obs Method  Exptime Qual Skysb Airmass Reducer\n")

    objects = []

    for f in flist:
        objname = f.split('_')[1].split('.')[0]

        if 'PTF' not in objname:
            continue
        if '_A_' in f or '_B_' in f:
            continue

        sp = np.load(f)[0]
        if 'quality' in sp:
            qual = sp['quality']
            if qual >= 3:
                continue
            else:
                objects.append(objname.replace("PTF", ""))
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
            if len(f.split('_')) > 2:
                obs = f.split('_')[-1].split('.')[0]
            elif len(f.split('_')) == 2:
                obs = f.split('_')[-1].split('.')[0]
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
        # Don't count missing objects
        if qual < 4:
            totexpt += expt
        else:
            lostexp += expt

        out.write("%-25s %4s %6s   %6.1f %4d %5s  %5.3f   %s\n" %
                  (objname, obs, meth, expt, qual, ("on" if skysub else "off"),
                   air, reducer))

    if len(objects) > 0:
        out.write("\n")
        out.write("Spectra are available in the marshal: \n")
        for o in objects:
            out.write("http://ptf.caltech.edu/cgi-bin/ptf/transient/view_source.cgi?name=%s\n"%(o))
        out.close()
        return True
    else:
        out.close()
        return False

if __name__ == '__main__':
    report = report()
    if report:
        current_dir = os.path.basename(os.path.abspath("."))
        # changed to iptftransient@lists.astro.caltech.edu on 9/24/2016
        cmd = 'cat report_ptf.txt | mail -s "SEDM DRP Report for %s" iptftransient@lists.astro.caltech.edu'%current_dir
        subprocess.call(cmd, shell=True)
    else:
        print("No PTF objects to report. \n")
