"""Divide a science frame by the flat frame with flexure

This routine may be obsolete and is currently unused by the pipeline.

Note:
    This is used as a python script as follows::

        usage: DivideFlat.py [-h] [--flexnpy FLEXNPY] [--outfile OUTFILE]
                     flatfits toflattenfits

        positional arguments:
          flatfits           Flat field fits file
          toflattenfits      To flatten fits file

        optional arguments:
          -h, --help         show this help message and exit
          --flexnpy FLEXNPY  Flexure .npy file
          --outfile OUTFILE  Output filename

"""

import argparse
import numpy as np
import astropy.io.fits as pf

import NPK.Fit as FF
import SEDMr.IO as IO

import SEDMr.Extraction as Extraction
import SEDMr.Wavelength as Wavelength
reload(FF)
reload(Extraction)
reload(Wavelength)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Divides a science frame by the flat frame, handles flexure

        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('flatfits', type=str, help='Flat field fits file')
    parser.add_argument('toflattenfits', type=str, help='To flatten fits file')
    parser.add_argument('--flexnpy', type=str, help='Flexure .npy file', default=None)

    parser.add_argument('--outfile', type=str, help='Output filename', default=None)

    args = parser.parse_args()

    flat = pf.open(args.flatfits)[0].data
    toflat = pf.open(args.toflattenfits)
    if args.flexnpy is not None:
        flex = np.load(args.flexnpy)[0]
        nmToPix = flex['skyline']/240
        dX = np.int(np.round(flex['dXnm'] * nmToPix))
        dY = np.int(np.round(flex['dYpix']))
        print "Rolling by %i/%i" % (dX, dY)
        flat = np.roll(flat, dY, 0)
        flat = np.roll(flat, dX, 1)
        toflat[0].header["FLATROLL"] = ("%s/%s" % (dX, dY), 
                                        "Flat field flexure correction")

    toflat[0].data /= flat
    toflat[0].header["FLATBY"] = (args.toflattenfits, "Flat field applied")
    toflat[0].data = toflat[0].data.astype(np.float64)
    IO.writefits( toflat, args.outfile, clobber=True, no_lossy_compress=True)

