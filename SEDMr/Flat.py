
import argparse
import numpy as np
import pylab as pl
import warnings

from SEDMr.FlatCorrection import FlatCorrection as FC


def measure_flat(extraction, fmeta, lamstart=700, lamend=900):

    corrections = []
    xs = []
    ys = []
    for i, e in enumerate(extraction):
        fc = FC(seg_id=e.seg_id)
        corrections.append(fc)

        if not e.ok:
            continue

        try:
            l, f = e.get_flambda()
        except:
            continue

        xs.append(np.nanmin(e.xrange)+e.xrefpix)
        ys.append(np.mean(e.yrange))

        roi = (l > lamstart) & (l <= lamend)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            fc.correction = np.nanmean(f[roi])

    vals = [f.get_correction(0) for f in corrections]
    medval = np.nanmedian(vals)

    ss = []
    for c in corrections: 
        try: 
            c.correction /= medval
            if c.correction > 2:
                c.correction = 1.0
            ss.append(c.correction)
        except:
            pass
    
    pl.figure(1)
    pl.clf()
    pl.scatter(xs, ys, c=ss, s=35, linewidth=0, vmin=0.8, vmax=1.2, marker='h',
               cmap=pl.get_cmap('jet'))
    pl.xlim(-100, 2048+100)
    pl.ylim(-100, 2048+100)
    pl.colorbar()
    pl.xlabel("X pixel @ %6.1f nm" % fmeta['fiducial_wavelength'])
    pl.ylabel("Y pixel")
    pl.title("Correction from %s to %s nm from %s" % (lamstart, lamend,
                                                      fmeta['outname']))
    if 'drp_version' in fmeta:
        drp_ver = fmeta['drp_version']
        ax = pl.gca()
        ax.annotate('DRP: ' + drp_ver, xy=(0.0, 0.01), xytext=(0, 0),
                    xycoords=('axes fraction', 'figure fraction'),
                    textcoords='offset points', size=10,
                    ha='center', va='bottom')
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
