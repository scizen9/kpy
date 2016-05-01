"""Graphical check of either a spectrum or a cube.

Generate matplotlib plots for visual confirmation of
extractions, or generate PDF files of the extractions.
Also can produce an ascii spectrum suitable for upload
or input to transient classification programs.

Functions
    :func:`checkSpec`   Plot a spectrum
    :func:`checkCube`   Plot a cube

"""
import argparse
import datetime
import os
import sys

import numpy as np
import pylab as pl
import scipy.signal
from scipy.stats import sigmaclip

import IO
import NPK.Standards as SS


def checkSpec(specname, corrname='std-correction.npy',
        redshift=0, smoothing=0, savefig=False, savespec=False):
    """Plot a spectrum.

    This function can produce an ascii spectrum suitable for uploading
    to the iPTF marshal and including header lines that describe the
    spectrum.

    Args:
        specname (str): name of spectrum numpy file, typically prefixed with "sp_"
        corrname (str): name of the calibration correction numpy file
        redshift (float): readshift for wavelength scale
        smoothing (int): number of pixels to smooth over
        savefig (bool): save a pdf of the plot
        savespec (bool): save an ascii spectrum

    Returns:
        None

    """

    if not os.path.isfile(specname):
        sys.exit("No such file: %s" % specname)

    # IO.readspec applies the calibration in the file specified
    lam, spec, skyspec, stdspec, ss, meta = \
            IO.readspec(specname,corrname=corrname)
    
    # Convert to Angstrom sized bins
    lam *= 10.
    spec /= 10.
    skyspec /= 10.
    stdspec /= 10.

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
        date = datetime.datetime(int(parts[0]), 1, 1) + \
                datetime.timedelta(int(parts[1])-1)
        utc = date.strftime("%Y %m %d") + " %s:%s:%s" % \
                (parts[2],parts[3],parts[4])
    except: utc = ''

    try: obj = hdr['OBJECT'].split()[0]
    except: obj = ''

    # Annotate plots
    pl.title("%s\n(airmass: %1.2f | Exptime: %i)" % 
            (specname, ec, et))
    pl.xlabel("Wavelength [Ang]")
    pl.ylabel("erg/s/cm2/ang")

    # Handle plot geometry
    plm = pl.get_current_fig_manager()
    plm.window.wm_geometry("900x500+10+10")

    # See if this is a standard star
    pred = specname.lstrip("sp_STD-")
    pred = pred.split("_")[0]
    pred = pred.lower().replace("+","").replace("-","_")
    if pred in SS.Standards:

        # Remove nm to Ang conversion (?)
        spec *= 10.
        skyspec *= 10.
        stdspec *= 10.

        # Get reference spectrum
        standard = SS.Standards[pred]
        slam = standard[:,0]
        sflx = standard[:,1] * 1.e-16

        # Calculate ratio in select region of spectrum
        lroi = (lam > 4500) & (lam < 6500)
        lmed = np.median(spec[lroi])

        sroi = (slam > 4500) & (slam < 6500)
        smed = np.median(sflx[sroi])

        rat = (lmed/smed)

        # Report offset
        ratmag = 2.5*np.log10(rat)
        print "Ref offset (Ref - Obs): %6.2f mag" % ratmag

        # Apply offset for plotting
        spec = spec / rat

    # Set wavelength range
    OK = (lam > 3800) & (lam < maxwl)

    # Apply redshift
    lamz = lam/(1+redshift)

    # Plot object spectrum
    
    # Plot limits
    pl.xlim(3600,maxwl+200)

    mx = np.nanmax(spec[OK])
    pl.ylim(-mx/10,mx+(mx/20))

    # No smoothing
    if smoothing == 0:
        pl.step(lamz[OK], spec[OK], linewidth=3)    # Plot data
    # Smoothing
    else:
        if smoothing > 5: order = 2
        else: order = 1
        smoothed = scipy.signal.savgol_filter(spec[OK], smoothing, order)
        pl.step(lamz[OK], smoothed, linewidth=3)    # Plot smoothed data

    # Legend for plot
    legend = ["obj",]

    # Overplot sky spectrum
    if skyspec is not None:
        pl.step(lamz[OK], skyspec[OK])
        legend.append("sky")

    # Overplot standard deviation spectrum
    if stdspec is not None:
        pl.step(lamz[OK], stdspec[OK])
        legend.append("err")

    # Overplot reference spectrum
    if pred in SS.Standards:
        print "Overplotting %s reference spectrum" % pred
        legend.append("ref")
        pl.plot(slam, sflx)

        # Remove plotting offset
        spec = spec * rat

    # Add legend
    pl.legend(legend)

    pl.grid(True)
    pl.ioff()

    # Save fig to file
    if savefig:
        outf = specname[(specname.find('_')+1):specname.find('.')]+'_SEDM.pdf'
        pl.savefig(outf)
        print "Figure saved to "+outf
    else:
        pl.show()

    if savespec:
        roi = (lam > 3800.0) & (lam < maxwl)
        wl = lam[roi]
        fl = spec[roi]
        srt = wl.argsort().argsort()
        outf = specname[(specname.find('_')+1):specname.find('.')]+'_SEDM.txt'
        header = "TELESCOPE: P60\nINSTRUMENT: SED-Machine\nUSER: %s" % user
        header = header + "\nOBJECT: %s\nOUTFILE: %s" % (obj, outf)
        header = header + "\nOBSUTC: %s\nEXPTIME %i" % (utc, et)
        header = header + "\nAIRMASS: %1.2f" % ec
        np.savetxt(outf, np.array([wl[srt], fl[srt]]).T, fmt='%8.1f  %.4e', 
                    header=header)
        print "Saved to "+outf


def checkCube(cubename, showlamrms=False, savefig=False):
    """Plot a datacube.

        This function can produce a PDF file of the data cube.

        Args:
            cubename (str): name of cube numpy file
            showlamrms (bool): display the rms of the wavelength fit
            savefig (bool): save a pdf of the plot

        Returns:
            None

        """
    if not os.path.isfile(cubename):
        sys.exit("No such file: %s" % cubename)

    cc, meta = np.load(cubename)
    fid_wave = meta['fiducial_wavelength']

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
        print("Nspax: %d, Nclip: %d, <RMS>: %f, RMS(std): %f" % 
                    (len(cc), (len(cc)-len(c)), Smdn, Sstd))
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
    pl.xlabel("X position [as] @ %6.1f nm" % fid_wave)
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
        """Check.py

        """, formatter_class=argparse.RawTextHelpFormatter)


    parser.add_argument('--cube', type=str, help='Fine correction path')
    parser.add_argument('--lambdarms', action="store_true", default=False, 
            help='Show lambda soln rms')
    parser.add_argument('--savefig', action="store_true", default=False, 
            help='Save pdf figure')
    parser.add_argument('--savespec', action="store_true", default=False, 
            help='Save spec ASCII file')
    parser.add_argument('--spec', type=str, help='Extracted spectrum file')
    parser.add_argument('--corrname', type=str, default='std-correction.npy')
    parser.add_argument('--redshift', type=float, default=0, help='Redshift')
    parser.add_argument('--smoothing', type=float, default=0, 
            help='Smoothing in pixels')

    args = parser.parse_args()

    if args.cube is not None:
        checkCube(args.cube, showlamrms=args.lambdarms, savefig=args.savefig)
    if args.spec is not None:
        checkSpec(args.spec, corrname=args.corrname, redshift=args.redshift,
            smoothing=args.smoothing, 
            savefig=args.savefig, savespec=args.savespec)


