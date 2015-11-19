
import argparse
import numpy as np
import os
import pylab as pl
import pyfits as pf
from scipy.interpolate import interp1d
import scipy.signal
import sys

import IO

 
def checkSpec(specname, corrname='std-correction.npy', redshift=0, smoothing=0):

    if not os.path.isfile(corrname) and corrname is not None:
        print "Loading old standard correction"
        corrname = '/scr2/npk/sedm/OUTPUT/2015mar25/std-correction.npy'
    else:
	print "Using current standard correction"

    if corrname is not None:
        corr = np.load(corrname)[0]
        corf = interp1d(corr['nm'], corr['correction'], bounds_error=False,
            fill_value=0.0)
    else:
        corf = lambda x: 1.0
        
    corf = lambda x: 1.0


    lam, spec, skyspec, stdspec, ss, meta = IO.readspec(specname)
    lam = np.roll(lam, 3)

    print "Plotting spectrum in %s" % specname
    try: print "Extraction radius: %1.2f" % ss['radius_as']
    except: pass

    try: ec = meta['airmass']
    except: ec = 0

    try: et = ss['exptime']
    except: et = 0

    pl.title("%s\n(airmass: %1.2f | Exptime: %i)" % (specname, ec, et))
    pl.xlabel("Wavelength [nm]")
    pl.ylabel("erg/s/cm2/ang")

    OK = (lam > 380) & (lam < 1000)
    legend = ["obj",]
    lamz = lam/(1+redshift)
    if smoothing == 0:
        pl.step(lamz[OK], spec[OK]*corf(lam[OK]), linewidth=3)
    else:
        if smoothing > 5: order = 2
        else: order = 1
        smoothed = scipy.signal.savgol_filter(spec[OK], smoothing, order)
        pl.step(lamz[OK], smoothed*corf(lamz[OK]), linewidth=3)

    if skyspec is not None:
        pl.step(lamz[OK], skyspec[OK]*corf(lam[OK]))
        legend.append("sky")

    if stdspec is not None:
        pl.step(lamz[OK], stdspec[OK]*corf(lam[OK]))
        legend.append("err")

    pl.legend(legend)
    pl.xlim(360,1000)

    roi = (lam > 400) & (lam < 999)
    mx = np.max(spec[roi]*corf(lam[roi]))
    pl.ylim(-mx/10,mx)
    pl.grid(True)
    pl.ioff()
    pl.show()

    wl = lam[roi]*10.
    fl = spec[roi]*corf(lam[roi])
    srt = wl.argsort().argsort()
    outf = specname[(specname.find('_')+1):specname.find('.')]+'_SEDM.txt'
    np.savetxt(outf, np.array([wl[srt], fl[srt]]).T)
    print "saved to "+outf


def checkCube(cubename, showlamrms=False):
    ''' Plot a datacube for checking'''
    
    cc = np.load(cubename)

    Xs = [c.X_as for c in cc]
    Ys = [c.Y_as for c in cc]
    if showlamrms:
	Ss = [0.] * len(cc)
	for i in range(len(cc)):
	    if cc[i].lamrms is not None:
		Ss[i] = cc[i].lamrms * 100.

    	smx = 0.1
    	smn = 0.0
	cbtitle = "Wavelength RMS [%]"
    else:
    	Ss = [c.trace_sigma for c in cc]
	smx = 2
	smn = 0.8
	cbtitle = "RMS trace width [pix]"

    pl.figure(1)
    pl.scatter(Xs, Ys, marker='H', linewidth=0, s=50, c=Ss, vmin=smn, vmax=smx)
    pl.title("Hexagonal Grid of Cube Positions")
    pl.xlim(-25,25)
    pl.ylim(-25,25)

    pl.colorbar(label=cbtitle)
    pl.xlabel("X position [as]")
    pl.ylabel("Y position [as]")
    pl.grid(True)
    pl.ioff()
    pl.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''Check.py

        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('--cube', type=str, help='Fine correction path')
    parser.add_argument('--lambdarms', action="store_true", default=False, help='Show lambda soln rms')
    parser.add_argument('--spec', type=str, help='Extracted spectrum file')
    parser.add_argument('--corrname', type=str, default='std-correction.npy')
    parser.add_argument('--redshift', type=float, default=0, help='Redshift')
    parser.add_argument('--smoothing', type=float, default=0, help='Smoothing in pixels')


    args = parser.parse_args()

    if args.cube is not None:
        checkCube(args.cube, showlamrms=args.lambdarms)
    if args.spec is not None:
        checkSpec(args.spec, corrname=args.corrname, redshift=args.redshift,
            smoothing=args.smoothing)


