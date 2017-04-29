"""Subtract background from ifu images."""

import argparse
import os
import numpy as np
import astropy.io.fits as pf
import sys
import shutil

import SEDMr.IO as IO

from scipy.ndimage.filters import gaussian_filter


def estimate_background(fine, infile, gausswidth=100, outname=None,
                       outint=False, fft_filt=False):

    if outname is None:
        print("Need an output name")
        return

    infile[0].data = infile[0].data.astype(np.float64)
    data = infile[0].data.copy()

    # loop over each trace
    for ff in fine:
        if not ff.ok and not ff.bkg_ok:
            continue
        if ff.xrange is None:
            continue
        if ff.poly is None:
            continue

        # get trace spatial ranges (xs, ys)
        xs = np.arange(*ff.xrange)
        ys = np.round(np.poly1d(ff.poly)(xs)).astype(np.int)

        # mask above and below each trace with nan's
        for dY in range(-4, 5):
            ty = ys.copy() - dY
            try:
                data[ty, xs] = np.nan
            except:
                pass

    from astropy.convolution import convolve, convolve_fft, Box2DKernel

    print("Traditional convolve (pass 1)")
    # get convolution kernel
    k = Box2DKernel(17)
    # start with original masked image for background
    bkg = data.copy()
    # write out starting image if requested
    if outint:
        IO.writefits(pf.PrimaryHDU(bkg), "test_0.fits", clobber=True)
        print("Wrote test_0.fits.gz")
    # keep track of nans
    nans = ~np.isfinite(data)
    # good data
    oks = np.isfinite(data)
    # iterate five times to remove object light from bkg
    for i in range(5):
        # convolve entire image
        bkg = convolve(bkg, k)
        # replace background light
        bkg[oks] = data[oks]
        print("\tIteration %d of 5" % (i+1))
        # write out each iteration if requested
        if outint:
            IO.writefits(pf.PrimaryHDU(bkg), "test_%i.fits" % (i+1),
                         clobber=True)
            print("Wrote test_%i.fits.gz" % i)

    # insert iterative smoothed object pixels
    data[nans] = bkg[nans]

    # any remaining nans?
    nans = ~np.isfinite(data)
    # fill them in with median value
    if np.count_nonzero(nans) > 0:
        data[nans] = np.nanmedian(data[oks])
    # write out intermediate result if requested
    if outint:
        fname = os.path.join(os.path.dirname(outname),
                             "lf_" + os.path.basename(outname))
        IO.writefits(data, fname, clobber=True)
        print("Wrote %s" % fname + ".gz")

    # use FFT filter if requested
    if fft_filt:
        print("FFT convolve (pass 2)")
        k = Box2DKernel(70)
        bkg = convolve_fft(data, k)
    # else, use gaussian filter of requested width in pixels
    else:
        print("Gaussian filter with width = %d (pass 2)" % gausswidth)
        bkg = gaussian_filter(data, gausswidth)

    # write resulting background to a gzipped fits file
    ofname = os.path.join(os.path.dirname(outname),
                          "bgd_" + os.path.basename(outname))
    HDU = pf.PrimaryHDU(bkg)
    HDU.header["GAUFWID"] = (gausswidth, 'Gaussian filter width in pixels')
    IO.writefits(HDU, ofname, clobber=True)
    print("Background image in %s" % ofname + ".gz")

    # record which file in output header
    infile[0].header["BGDSUB"] = "Background subtracted using %s" % ofname
    infile[0].header["GAUFWID"] = (gausswidth, 
                                   'Gaussian filter width in pixels')
    ofname = os.path.join(os.path.dirname(outname),
                          "bs_" + os.path.basename(outname))
    # subtract background
    infile[0].data -= bkg
    # write out the resulting fits file
    IO.writefits(infile, ofname, clobber=True)
    print("Subtracted image in %s" % ofname + ".gz")

    return bkg


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        '''SubtractBackground.py

        ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('fine', type=str, help='Fine correction path')
    parser.add_argument('infile', type=str, help='Path to FITS file to refit')
    parser.add_argument('--gausswidth', type=int, default=100,
                        help='Gaussian filter width in pixels')
    parser.add_argument('--outint', default=False, action="store_true",
                        help='Output intermediate images')
    parser.add_argument('--fftfilt', default=False, action="store_true",
                        help='Use FFT convolve instead of Gaussian filter')

    args = parser.parse_args()
    fine = np.load(args.fine)
    infile = pf.open(args.infile)

    if infile[0].header['EXPTIME'] < 30:
        fname = os.path.join(os.path.dirname(args.infile),
                             "bs_" + os.path.basename(args.infile))
        shutil.copy(args.infile, fname)
        os.system("gzip --fast --force %s" % fname)
        sys.exit(0)

    gauss_width = args.gausswidth

    background = estimate_background(fine, infile, gausswidth=gauss_width,
                                     outname=args.infile, outint=args.outint,
                                     fft_filt=args.fftfilt)
