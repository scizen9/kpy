# -*- coding: utf-8 -*-
"""
Created on Thu Aug 11 11:23:00 2016

@author: nblago
"""
import datetime
import glob
import os
import argparse
import re
import SEDMr.RunSnid as RunSnid


def classify(spec_dir='./', overwrite=False):
    """
    Runs snid in batch mode on all the *_SEDM.txt files found in the given
    directory.  If a given file was already classified, it skips it, unless
    overwrite is requested.
    """

    summary = []
    flo = glob.glob(os.path.join(spec_dir, "*_SEDM.txt"))
    fln = glob.glob(os.path.join(spec_dir, "spec_*.txt"))
    flb = flo + fln
    for fl in flb:
        print(fl)
        # don't classify standard stars
        if "STD" in fl or "BD" in fl or "Feige" in fl or "HZ" in fl:
            print("standard star")
            continue
        if "Hiltner" in fl or "Kopff" in fl or "LTT" in fl:
            print("standard star")
            continue
        # don't classify galaxies
        if "NGC" in fl or "PGC" in fl or "MCG" in fl or "2MASX" in fl:
            print("galaxy")
            continue
        # don't classify stars
        if "TYC" in fl or "SAO" in fl or "HD" in fl or "Tycho" in fl:
            print("star")
            continue
        # skip uncalibrated objects
        if "notfluxcal" in fl:
            print("uncalibrated")
            continue
        # retrieve the quality of the spectra.
        with open(fl, "r") as sfl:
            lines = sfl.readlines()

            q = [li for li in lines if "QUALITY" in li]

            if len(q) > 0:
                token = re.search(r'\(?([0-9]+)\)?', q[0])
                q = int(token.group(1))
            else:
                if "crr_b_ifu" in fl:
                    q = 1
                else:
                    q = 5

            # If quality is good, check for previous classification
            if q < 3:
                clas = [li for li in lines if "SNID" in li]
                # If the file has been classified, move to the next
                if len(clas) > 0 and not overwrite:
                    print("already classified")
                    continue
            else:
                print("low quality")
                continue
            
        # If we are here, we run the classification with snid
        specfl = fl.split('/')[-1]
        good = RunSnid.run_snid(spec_file=specfl, overwrite=overwrite)
        # If we actually ran, record the results
        if good:
            res = specfl + " " + RunSnid.record_snid(spec_file=specfl)
            summary.append(res)
    # END loop over each file matching *_SEDM.txt
    for res in summary:
        print(res)

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        """

        Copies the output ascii spectra to the iPTF marshal

        Checks to be sure spectra have quality > 3 and only
        copies spectra for objects with "PTF" prefix.

        Specify input directory with -d, or --specdir parameters.
        If none specified, use current date directory in /scr2/sedmdrp/redux/

        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('-d', '--specdir', type=str, dest="specdir",
                        help='Directory with output PTF*_SEDM.txt files.',
                        default=None)
    parser.add_argument('--overwrite', action="store_true", default=False,
                        help='Overwrite existing classification')

    args = parser.parse_args()
    
    specdir = args.specdir
    
    if specdir is None:
        timestamp = datetime.datetime.isoformat(datetime.datetime.utcnow())
        timestamp = timestamp.split("T")[0].replace("-", "")
        specdir = os.path.join("/scr2/sedmdrp/redux/", timestamp)
    else:
        specdir = os.path.abspath(specdir)
        timestamp = os.path.basename(specdir)
        os.chdir(specdir)

    print(os.getcwd())

    # Run snid on extracted spectra
    classify(spec_dir=specdir, overwrite=args.overwrite)
