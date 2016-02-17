
import argparse
import pdb
import numpy as np
import pylab as pl
import pyfits as pf
import sys
import warnings


import NPK.Fit as FF
from astropy.table import Table 


from scipy.spatial import KDTree 
import scipy.signal as SG
from scipy.interpolate import interp1d


from numpy.polynomial.chebyshev import chebfit, chebval

import SEDMr.Extraction as Extraction
import SEDMr.Wavelength as Wavelength
from SEDMr.FlatCorrection import FlatCorrection as FC

reload(FF)
reload(Extraction)
reload(Wavelength)


def measure_flat(extraction, meta, 
        lamstart=700,
        lamend=900,
        fidwave=None,
        outfile='flat.npy'):

    if fidwave is None:
        fid_wave = Wavelength.fiducial_wavelength()
    else:
        fid_wave = fidwave


    corrections = []
    Xs = []
    Ys = []
    for i,e in enumerate(extraction):
        fc = FC(seg_id=e.seg_id)
        corrections.append(fc)

        if not e.ok: continue

        try: l,f = e.get_flambda()
        except: continue

        X = np.argmin(np.abs(l-fid_wave))
        Xs.append(e.xrange[0] + X)
        Ys.append(np.mean(e.yrange))

        ROI = (l>lamstart) & (l <= lamend)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            fc.correction = np.nanmean(f[ROI])

    vals = [f.get_correction(0) for f in corrections]
    medval = np.median(vals)

    Ss = []
    for c in corrections: 
        try: 
            c.correction /= medval
            if c.correction > 2:
                c.correction = 1.0
            Ss.append(c.correction)
        except: pass
    
    pl.figure(1)
    pl.clf()
    pl.scatter(Xs, Ys, c=Ss,s=70,linewidth=0, vmin=0.8,vmax=1.2,marker='h')
    pl.xlim(-100,2048+100)
    pl.ylim(-100,2048+100)
    pl.colorbar()
    pl.xlabel("X pixel @ %6.1f nm" % fid_wave)
    pl.ylabel("Y pixel")
    pl.title("Correction from %s to %s nm from %s" % (lamstart, lamend,
                meta['outname']))
    pl.savefig("flat-field-values.pdf")

    return corrections

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        """Create dome flat

        """, formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('infile', type=str, help='Path to dome flat')
    parser.add_argument('--lamstart', type=float,
            help='Wavelength range start', default=700.0)
    parser.add_argument('--lamend', type=float, 
            help='Wavelength range end', default=900.0)
    parser.add_argument('--outfile', type=str, 
            help='Output filename', default="flat-dome-700to900.npy")

    args = parser.parse_args()

    ext, meta = np.load(args.infile)
    if 'fiducial_wavelength' in meta:
        fidwave = meta['fiducial_wavelength']
    else:
        fidwave = None
    flat = measure_flat(ext, meta, lamstart=args.lamstart, lamend=args.lamend,
                        fidwave=fidwave)

    np.save(args.outfile, flat)
    print "Wrote %s" % args.outfile
