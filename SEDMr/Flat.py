
import argparse
import numpy as np
import pylab as pl
import warnings

import NPK.Fit as FF

import SEDMr.Extraction as Extraction
import SEDMr.Wavelength as Wavelength
from SEDMr.FlatCorrection import FlatCorrection as FC

reload(FF)
reload(Extraction)
reload(Wavelength)


def measure_flat(extraction, fmeta, lamstart=700, lamend=900):

    corrections = []
    Xs = []
    Ys = []
    for i, e in enumerate(extraction):
        fc = FC(seg_id=e.seg_id)
        corrections.append(fc)

        if not e.ok:
            continue

        try:
            l, f = e.get_flambda()
        except:
            continue

        Xs.append(np.nanmin(e.xrange)+e.xrefpix)
        Ys.append(np.mean(e.yrange))

        ROI = (l > lamstart) & (l <= lamend)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            fc.correction = np.nanmean(f[ROI])

    vals = [f.get_correction(0) for f in corrections]
    medval = np.nanmedian(vals)

    Ss = []
    for c in corrections: 
        try: 
            c.correction /= medval
            if c.correction > 2:
                c.correction = 1.0
            Ss.append(c.correction)
        except:
            pass
    
    pl.figure(1)
    pl.clf()
    pl.scatter(Xs, Ys, c=Ss, s=35, linewidth=0, vmin=0.8, vmax=1.2, marker='h',
               cmap=pl.get_cmap('jet'))
    pl.xlim(-100, 2048+100)
    pl.ylim(-100, 2048+100)
    pl.colorbar()
    pl.xlabel("X pixel @ %6.1f nm" % fmeta['fiducial_wavelength'])
    pl.ylabel("Y pixel")
    pl.title("Correction from %s to %s nm from %s" % (lamstart, lamend,
                                                      fmeta['outname']))
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
    flat = measure_flat(ext, meta, lamstart=args.lamstart, lamend=args.lamend)

    np.save(args.outfile, flat)
    print("Wrote %s" % args.outfile)
