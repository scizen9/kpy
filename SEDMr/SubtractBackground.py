"""Subtract background from ifu images."""

import argparse
import os
import numpy as np
import pyfits as pf
import sys
import shutil

import SEDMr.IO as IO

from scipy.ndimage.filters import gaussian_filter

from scipy.weave import converters
import scipy.weave as weave
 
 
def weaveConvolve(image, kernel):
 
    image = image.astype(float)
    kernel = kernel.astype(float)
 
    nx, ny = image.shape
    nkx, nky = kernel.shape
 
    if nkx % 2 == 0 or nky % 2 == 0:
        raise Exception("Kernel dimensions should be odd")
 
    smoothed = np.zeros(image.shape)
    isvalid = ~np.isnan(image)
 
    code = """
            double top, bot;
            int wkx, wky, iimin, iimax, jjmin, jjmax;
            wkx = (nkx-1)/2;
            wky = (nky-1)/2;
            for (int i=0; i<nx; ++i) {
                for (int j=0; j<ny; ++j) {
                    if(isvalid(i,j)) {
                        top = 0.;
                        bot = 0.;
                        if(i-wkx > 0) { iimin = i-wkx; } else { iimin = 0; };
                        if(i+wkx < nx-1) { iimax = i+wkx; } else { iimax = nx-1; };
                        if(j-wkx > 0) { jjmin = j-wky; } else { jjmin = 0; };
                        if(j+wkx < ny-1) { jjmax = j+wky; } else { jjmax = ny-1; };
                        for (int ii=iimin; ii <= iimax ; ++ii) {
                            for (int jj=jjmin; jj <= jjmax; ++jj) {
                                if(isvalid(ii,jj)) {
                                    top = top + kernel(wkx + ii-i,wky + jj-j) * image(ii,jj);
                                    bot = bot + kernel(wkx + ii-i,wky + jj-j);
                                }
                            }
                        }
                        smoothed(i,j) = top / bot;
                    } else {
                        smoothed(i,j) = image(i,j);
                    }
                }
            }
            return_val = 1;
            """
 
    weave.inline(code, ['image', 'isvalid', 'nx', 'ny', 'kernel', 'nkx', 'nky', 'smoothed'],
                 type_converters=converters.blitz, compiler = 'gcc')
 
    return smoothed


def remove(fname):
    try: os.remove(fname)
    except: pass



def estimateBackground(fine, infile, flex=None, gausswidth=100, outname=None):

    if outname is None:
        print "Need an output name"
        return

    dX = 0
    dY = 0

    infile[0].data = infile[0].data.astype(np.float64)
    data = infile[0].data.copy()
    

    for ff in fine:
        if not ff.ok: continue
        if ff.xrange is None: continue
        if ff.poly is None: continue
        
        xs = np.arange(*ff.xrange)
        ys = np.round(np.poly1d(ff.poly)(xs)).astype(np.int)

        for dY in xrange(-5,6):
            ty = ys.copy() - dY
            try: data[ty,xs] = np.nan
            except: pass

    from astropy.convolution import convolve, convolve_fft, Box2DKernel

    print "Traditional convolve (pass 1)"
    k = Box2DKernel(17)
    flt = data.copy()
    NaNs = ~np.isfinite(data)
    OKs = np.isfinite(data)
    for i in xrange(5):
        flt = convolve(flt, k)
        flt[OKs] = data[OKs]
        print "\tIteration %d of 5" %  (i+1)
        #IO.writefits(pf.PrimaryHDU(flt), "test_%i.fits.gz" % i, clobber=True)

    data[NaNs] = flt[NaNs]
    fname = os.path.join(os.path.dirname(outname), 
        "lf_" + os.path.basename(outname))
    #IO.writefits(data, fname, clobber=True)
    

    #print "FFT convolve (pass 2)"
    print "Gaussian filter with width = %d (pass 2)" % gausswidth
    #k = Box2DKernel(70)
    #flt = convolve_fft(data, k)
    flt = gaussian_filter(data, gausswidth)


    fname = os.path.join(os.path.dirname(outname), 
        "bgd_" + os.path.basename(outname))
    HDU = pf.PrimaryHDU(flt)
    HDU.header["GAUFWID"] = (gausswidth, 'Gaussian filter width in pixels')
    IO.writefits(HDU, fname, clobber=True)
    print "Background image in %s" % fname + ".gz"
    

    infile[0].header["BGDSUB"] = "Background subtracted using %s" % fname
    infile[0].header["GAUFWID"] = (gausswidth, 
                                   'Gaussian filter width in pixels')
    fname = os.path.join(os.path.dirname(outname), 
        "bs_" + os.path.basename(outname))
    infile[0].data -= flt

    IO.writefits(infile, fname, clobber=True)
    print "Subtracted image in %s" % fname + ".gz"

    return flt



if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''SubtractBackground.py

        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('fine', type=str, help='Fine correction path')
    parser.add_argument('infile', type=str, help='Path to FITS file to refit')
    parser.add_argument('--flexfile', type=str, help='Path to flexure npy file')
    parser.add_argument('--gausswidth', type=int, default=100, help='Gaussian filter width in pixels')


    args = parser.parse_args()
    fine = np.load(args.fine)
    infile = pf.open(args.infile)

    if infile[0].header['EXPTIME'] < 30:
        fname = os.path.join(os.path.dirname(args.infile), 
            "bs_" + os.path.basename(args.infile))
        shutil.copy(args.infile, fname)
        os.system("gzip --fast --force %s" % fname)
        sys.exit(0)

    if args.flexfile is not None:
        flex = np.load(args.flexfile)
    else: flex=None

    gausswidth = args.gausswidth


    background = estimateBackground(fine, infile, flex=flex, 
                                    gausswidth=gausswidth, outname=args.infile)
