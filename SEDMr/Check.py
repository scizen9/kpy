
import argparse
import numpy as np
import os
import pylab as pl
import pyfits as pf
from scipy.interpolate import interp1d
import scipy.signal
from scipy.stats import sigmaclip
import sys
import datetime

import IO
import NPK.Standards as SS

 
def checkSpec(specname, redshift=0, smoothing=0, savefig=False):

    # IO.readspec applies the calibration in std-correction.npy
    lam, spec, skyspec, stdspec, ss, meta = IO.readspec(specname)
    
    lam = lam * 10.     # convert to Angstroms

    # Get header values
    hdr = meta['header']

    print "Plotting spectrum in %s" % specname
    try: print "Extraction radius: %1.2f asec" % ss['radius_as']
    except: pass

    try: ec = meta['airmass']
    except: ec = 0

    try: et = ss['exptime']
    except: et = 0

    try: maxwl = meta['maxnm'] * 10.
    except: maxwl = 9200.0

    print "Max Angstroms: %7.1f" % maxwl

    try: 
        ec2 = meta['airmass2']
        et *= 2
    except: ec2 = 0

    try: user = meta['user']
    except: user = ''

    try:
        utc = meta['utc']
        parts = utc.split(":")
        date = datetime.datetime(int(parts[0]), 1, 1) + datetime.timedelta(int(parts[1])-1)
        utc = date.strftime("%Y %m %d") + " %s:%s:%s" % (parts[2],parts[3],parts[4])
    except: utc = ''

    try: obj = hdr['OBJECT'].split()[0]
    except: obj = ''

    # Annotate plots
    pl.title("%s\n(airmass: %1.2f | Exptime: %i)" % (specname, ec, et))
    pl.xlabel("Wavelength [Ang]")
    pl.ylabel("erg/s/cm2/ang")

    # Handle plot geometry
    plm = pl.get_current_fig_manager()
    plm.window.wm_geometry("900x500+10+10")

    # Set plot limits
    OK = (lam > 3800) & (lam < maxwl)
    legend = ["obj",]
    lamz = lam/(1+redshift)

    # Smoothing option
    if smoothing == 0:
        pl.step(lamz[OK], spec[OK], linewidth=3)    # Plot data
    else:
        if smoothing > 5: order = 2
        else: order = 1
        smoothed = scipy.signal.savgol_filter(spec[OK], smoothing, order)
        pl.step(lamz[OK], smoothed, linewidth=3)    # Plot smoothed data

    # Overplot sky spectrum
    if skyspec is not None:
        pl.step(lamz[OK], skyspec[OK])
        legend.append("sky")

    # Overplot standard deviation spectrum
    if stdspec is not None:
        pl.step(lamz[OK], stdspec[OK])
        legend.append("err")

    # See if this is a standard star
    pred = specname.lstrip("sp_STD-")
    pred = pred.split("_")[0]
    pred = pred.lower().replace("+","").replace("-","_")
    if pred in SS.Standards:
        print "Overplotting %s reference spectrum" % pred
        legend.append("ref")

        standard = SS.Standards[pred]
        slam = standard[:,0]
        sflx = standard[:,1] * 1.e-16

        lroi = (lam > 4500) & (lam < 6500)
        lmed = np.median(spec[lroi])

        sroi = (slam > 4500) & (slam < 6500)
        smed = np.median(sflx[sroi])

        rat = (lmed/smed)
        ratmag = 2.5*np.log10(rat)
        print "Ref offset: %6.2f mag" % ratmag

        sflx = sflx * rat

        pl.plot(slam, sflx)

    pl.legend(legend)

    # Plot limits
    pl.xlim(3600,maxwl)

    roi = (lam > 3800.0) & (lam < maxwl)
    mx = np.nanmax(spec[roi])
    pl.ylim(-mx/10,mx)

    pl.grid(True)
    pl.ioff()

    # Save fig to file
    if savefig:
        outf = specname[(specname.find('_')+1):specname.find('.')]+'_SEDM.pdf'
        pl.savefig(outf)
        print "Figure saved to "+outf
    else:
        pl.show()

    wl = lam[roi]
    fl = spec[roi]
    srt = wl.argsort().argsort()
    outf = specname[(specname.find('_')+1):specname.find('.')]+'_SEDM.txt'
    header = 'TELESCOPE: P60\nINSTRUMENT: SED-Machine\nUSER: %s\nOBJECT: %s\nOUTFILE: %s\nOBSUTC: %s\nEXPTIME %i\nAIRMASS: %1.2f' % (user, obj, outf, utc, et, ec)
    np.savetxt(outf, np.array([wl[srt], fl[srt]]).T, fmt='%8.1f  %.4e', header=header)
    print "Saved to "+outf


def checkCube(cubename, showlamrms=False, savefig=False):
    ''' Plot a datacube for checking'''
    
    cc = np.load(cubename)

    Xs = [c.X_as for c in cc]
    Ys = [c.Y_as for c in cc]
    if showlamrms:
        Ss = [0.] * len(cc)
        for i in range(len(cc)):
            if cc[i].lamrms is not None:
                if np.isfinite(cc[i].lamrms):
                    Ss[i] = cc[i].lamrms

        c, low, upp = sigmaclip(Ss)
        Smdn = np.median(c)
        Sstd = np.nanstd(c)
        print "Nspax: %d, Nclip: %d, <RMS>: %f, RMS(std): %f" % (len(cc), (len(cc)-len(c)), Smdn, Sstd)
        smx = Smdn + 3.* Sstd
        smn = Smdn - 3.* Sstd
        if smn < 0.: smn = 0.
        cbtitle = "Wavelength RMS [nm]"
        outf = "cube_lambdarms.pdf"
    else:
        Ss = [c.trace_sigma for c in cc]
        smx = 2
        smn = 0.8
        cbtitle = "RMS trace width [pix]"
        outf = "cube_trace_sigma.pdf"

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
    if savefig:
        pl.savefig(outf)
        print "Figure saved to "+outf
    else:
        pl.show()

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''Check.py

        ''', formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('--cube', type=str, help='Fine correction path')
    parser.add_argument('--lambdarms', action="store_true", default=False, help='Show lambda soln rms')
    parser.add_argument('--savefig', action="store_true", default=False, help='Save pdf figure')
    parser.add_argument('--spec', type=str, help='Extracted spectrum file')
    parser.add_argument('--redshift', type=float, default=0, help='Redshift')
    parser.add_argument('--smoothing', type=float, default=0, help='Smoothing in pixels')


    args = parser.parse_args()

    if args.cube is not None:
        checkCube(args.cube, showlamrms=args.lambdarms, savefig=args.savefig)
    if args.spec is not None:
        checkSpec(args.spec, redshift=args.redshift,
            smoothing=args.smoothing, savefig=args.savefig)


