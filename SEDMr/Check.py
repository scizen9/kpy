"""Graphical check of either a spectrum or a cube.

Generate matplotlib plots for visual confirmation of
extractions, or generate PDF files of the extractions.
Also can produce an ascii spectrum suitable for upload
or input to transient classification programs.

Functions
    * :func:`check_cube`   Plot a cube
    * :func:`check_spec`   Plot a spectrum

Note:
    This is used as a python script as follows::

        usage: Check.py [-h] [--cube CUBE] [--lambdarms] [--savefig]
                [--savespec] [--spec SPEC] [--corrname CORRNAME]
                [--redshift REDSHIFT] [--smoothing SMOOTHING]

        optional arguments:
          -h, --help            show this help message and exit
          --cube CUBE           Fine correction path
          --lambdarms           Show lambda soln rms
          --savefig             Save pdf figure
          --savespec            Save spec ASCII file
          --spec SPEC           Extracted spectrum file
          --corrname CORRNAME
          --redshift REDSHIFT   Redshift
          --smoothing SMOOTHING
                                Smoothing in pixels
"""

import argparse
import datetime
import os
import sys
from builtins import input

import numpy as np
import pylab as pl
import scipy.signal
from scipy.stats import sigmaclip

import IO
import NPK.Standards as Stds


def check_cube(cubename, showlamrms=False, savefig=False, no_stamp=False):
    """Plot a datacube.

        This function can produce a PDF file of the data cube.

        Args:
            cubename (str): name of cube numpy file
            showlamrms (bool): display the rms of the wavelength fit,
                otherwise display average trace width
            savefig (bool): save a pdf of the plot
            no_stamp (bool): set to prevent printing DRP version stamp on plot

        Returns:
            None

        """

    if not os.path.isfile(cubename):
        sys.exit("No such file: %s" % cubename)

    cc, meta = np.load(cubename)
    fid_wave = meta['fiducial_wavelength']
    if 'drp_version' in meta:
        drp_ver = meta['drp_version']
    else:
        drp_ver = ''

    xs = [c.X_as for c in cc]
    ys = [c.Y_as for c in cc]
    if showlamrms:
        ss = [0.] * len(cc)
        for i in range(len(cc)):
            if cc[i].lamrms is not None:
                if np.isfinite(cc[i].lamrms):
                    ss[i] = cc[i].lamrms

        c, low, upp = sigmaclip(ss)
        smdn = np.nanmedian(c)
        sstd = np.nanstd(c)
        print("Nspax: %d, Nclip: %d, <RMS>: %f, RMS(std): %f" %
              (len(cc), (len(cc) - len(c)), float(smdn), float(sstd)))
        # smx = smdn + 3. * sstd
        # smn = smdn - 3. * sstd
        # if smn < 0.:
        #     smn = 0.
        smn = 0.
        smx = 1.2
        cbtitle = "Wavelength RMS [nm]"
        outf = "cube_lambdarms.pdf"
    else:
        ss = [c.trace_sigma for c in cc]
        smx = 2
        smn = 0.8
        cbtitle = "RMS trace width [pix]"
        outf = "cube_trace_sigma.pdf"

    pl.figure(1)
    pl.scatter(xs, ys, marker='H', linewidth=0, s=35, c=ss, vmin=smn, vmax=smx,
               cmap=pl.get_cmap('jet'))
    if not no_stamp:
        pl.title("Hexagonal Grid of Cube Positions")
    pl.xlim(15, -15)
    pl.ylim(-15, 15)

    pl.colorbar(label=cbtitle)
    pl.xlabel("RA offset [asec] @ %6.1f nm" % fid_wave)
    pl.ylabel("Dec offset [asec]")
    # Add drp version
    if len(drp_ver) > 0 and not no_stamp:
        ax = pl.gca()
        ax.annotate('DRP: '+drp_ver, xy=(0.0, 0.01), xytext=(0, 0),
                    xycoords=('axes fraction', 'figure fraction'),
                    textcoords='offset points', size=6,
                    ha='center', va='bottom')

    pl.grid(True)
    pl.ioff()
    if savefig:
        pl.savefig(outf)
        print("Figure saved to " + outf)
    else:
        pl.show()


def check_spec(specname, corrname='std-correction.npy', redshift=0, smoothing=0,
               savefig=False, savespec=False, interact=False, no_stamp=False):
    """Plot a spectrum.

    This function can produce an ascii spectrum suitable for uploading
    to the iPTF marshal and including header lines that describe the
    spectrum.

    Args:
        specname (str): name of spectrum numpy file, usually prefixed with "sp_"
        corrname (str): name of the calibration correction numpy file
        redshift (float): redshift for wavelength scale
        smoothing (int): number of pixels to smooth over
        savefig (bool): save a pdf of the plot
        savespec (bool): save an ascii spectrum
        interact (bool): query for the quality of the spectrum
        no_stamp (bool): set to prevent printing DRP version stamp on plot

    Returns:
        None

    """

    if not os.path.isfile(specname):
        sys.exit("No such file: %s" % specname)

    # IO.readspec applies the calibration in the file specified
    lam, spec, skyspec, stdspec, ss, meta = \
        IO.readspec(specname, corrname=corrname)

    # Convert to Angstrom sized bins
    lam *= 10.
    spec /= 10.
    if skyspec is not None:
        skyspec /= 10.
    if stdspec is not None:
        stdspec /= 10.

    # Print wavelength range
    print("Wavelengths from %.1f - %.1f" %
          (float(np.nanmin(lam[np.isfinite(spec)])),
           float(np.nanmax(lam[np.isfinite(spec)]))))

    # Get object name
    if 'header' in meta:
        hdr = meta['header']
        if 'OBJECT' in hdr:
            obj = hdr['OBJECT'].split()[0]
        else:
            obj = ''
    else:
        obj = ''

    print("Plotting spectrum in %s" % specname)
    if 'radius_as' in ss:
        print("Extraction radius: %1.2f asec" % ss['radius_as'])

    if 'airmass1' in meta:
        ec = meta['airmass1']
        if 'airmass2' in meta:
            ec = (ec + meta['airmass2']) / 2.
    elif 'airmass' in meta:
        ec = meta['airmass']
    else:
        ec = 0

    if 'exptime' in ss:
        et = ss['exptime']
    else:
        et = 0

    if 'maxnm' in meta:
        maxwl = meta['maxnm'] * 10.
    else:
        maxwl = 9200.0

    print("Max Angstroms: %7.1f" % maxwl)

    # If it has airmass2 then must be A/B pair
    if 'airmass2' in meta:
        et *= 2

    if 'user' in meta:
        user = meta['user']
    else:
        user = ''

    if 'sky_subtraction' in ss:
        skysub = ss['sky_subtraction']
    else:
        skysub = True

    if 'quality' in ss:
        qual = ss['quality']
    else:
        qual = 0

    try:
        utc = meta['utc']
        parts = utc.split(":")
        date = datetime.datetime(int(parts[0]), 1, 1) + \
               datetime.timedelta(int(parts[1]) - 1)
        utc = date.strftime("%Y %m %d") + " %s:%s:%s" % \
                                          (parts[2], parts[3], parts[4])
    except:
        utc = ''

    if 'drp_version' in ss:
        drp_ver = ss['drp_version']
    else:
        drp_ver = ''

    # Annotate plots
    if qual > 0:
        tlab = "%s\n(Air: %1.2f | Expt: %i | Skysub: %s | Qual: %d)" % \
             (specname, ec, et, "On" if skysub else "Off", qual)
    else:
        tlab = "%s\n(Air: %1.2f | Expt: %i | Skysub: %s)" % \
             (specname, ec, et, "On" if skysub else "Off")
    if not no_stamp:
        pl.title(tlab)
    pl.xlabel("Wavelength [Ang]")
    pl.ylabel("erg/s/cm2/ang")

    # Handle plot geometry
    # plm = pl.get_current_fig_manager()
    # plm.window.wm_geometry("900x500+10+10")

    # See if this is a standard star
    pred = specname[7:]
    pred = pred.split("_")[0]
    pred = pred.lower().replace("+", "").replace("-", "_")
    if pred in Stds.Standards:
        # Remove nm to Ang conversion (?)
        spec *= 10.
        skyspec *= 10.
        stdspec *= 10.

        # Get reference spectrum
        standard = Stds.Standards[pred]
        slam = standard[:, 0]
        sflx = standard[:, 1] * 1.e-16  # type: np.ndarray

        # Calculate ratio in select region of spectrum
        lroi = (lam > 4500) & (lam < 6500)
        lmed = np.nanmedian(spec[lroi])

        sroi = (slam > 4500) & (slam < 6500)
        smed = np.nanmedian(sflx[sroi])

        rat = (lmed / smed)

        # Report offset
        ratmag = 2.5 * np.log10(rat)
        print("Ref offset (Ref - Obs): %6.2f mag" % ratmag)

        # Apply offset for plotting
        spec /= rat

    # Set wavelength range
    ok = (lam > 3800) & (lam < maxwl)

    # Apply redshift
    lamz = lam / (1 + redshift)  # type: np.ndarray

    # Plot object spectrum

    # Plot limits
    pl.xlim(3600, maxwl + 200)

    mx = np.nanmax(spec[ok])
    pl.ylim(-mx / 10, mx + (mx / 20))

    # No smoothing
    if smoothing == 0:
        pl.step(lamz[ok], spec[ok], linewidth=3)  # Plot data
    # Smoothing
    else:
        if smoothing > 5:
            order = 2
        else:
            order = 1
        smoothed = scipy.signal.savgol_filter(spec[ok], smoothing, order)
        pl.step(lamz[ok], smoothed, linewidth=3)  # Plot smoothed data

    # Legend for plot
    legend = ["obj", ]

    # Overplot sky spectrum
    if skyspec is not None:
        pl.step(lamz[ok], skyspec[ok])
        legend.append("sky")

    # Overplot standard deviation spectrum
    if stdspec is not None:
        pl.step(lamz[ok], stdspec[ok])
        legend.append("err")

    # Overplot reference spectrum
    if pred in Stds.Standards:
        print("Overplotting %s reference spectrum" % pred)
        legend.append("ref")
        pl.plot(slam, sflx)

        # Remove plotting offset
        spec *= rat

    # Add legend
    pl.legend(legend)

    # Add drp version
    if len(drp_ver) > 0 and not no_stamp:
        ax = pl.gca()
        ax.annotate('DRP: '+drp_ver, xy=(0.0, 0.01), xytext=(0, 0),
                    xycoords=('axes fraction', 'figure fraction'),
                    textcoords='offset points', size=6,
                    ha='center', va='bottom')

    pl.grid(True)
    pl.ioff()

    # Get the culprit
    reducer = os.getenv("SEDM_USER")
    if reducer is None:
        reducer = "Unknown"
    # Load spectrum
    res = np.load(specname)

    if interact:
        # Set up plot
        pl.ion()
        pl.show()
        # Get quality of observation
        print("Enter quality of observation:")
        print("1 - good       (no problems)")
        print("2 - acceptable (minor problem)")
        print("3 - poor       (major problem)")
        print("4 - no object visible")
        q = 'x'
        qual = -1
        prom = ": "
        while qual < 1 or qual > 4:
            q = input(prom)
            if type(q) != int:
                if q.isdigit():
                    qual = int(q)
            else:
                qual = q
            if qual < 1 or qual > 4:
                prom = "Try again: "
        print("Quality = %d" % qual)
        tlab = "%s\n(Air: %1.2f | Expt: %i | Skysub: %s | Qual: %d)" % \
               (specname, ec, et, "On" if skysub else "Off", qual)
        pl.title(tlab)
        # Update quality
        res[0]['quality'] = qual

        # Get reducer
        print("Enter reducer of observations:")
        prom = "<cr> = ("+reducer+"): "
        q = input(prom)
        if len(q) > 0:
            reducer = q

        pl.ioff()

    # Add reducer and save spectrum
    res[0]['reducer'] = reducer
    print("Reducer: %s" % reducer)
    try:
        np.save(specname, res)
    except IOError:
        print("Unable to update %s with reducer" % specname)

    # Save fig to file
    if savefig:
        outf = specname[(specname.find('_') + 1):specname.rfind('.')] + \
               '_SEDM.pdf'
        pl.savefig(outf)
        print("Figure saved to " + outf)

    if not interact and not savefig:
        pl.show()

    if savespec:
        roi = (lam > 3800.0) & (lam < maxwl) & np.isfinite(spec)
        wl = lam[roi]
        fl = spec[roi]
        srt = wl.argsort().argsort()
        outf = specname[(specname.find('_') + 1):specname.rfind('.')] + \
            '_SEDM.txt'
        header = "TELESCOPE: P60\nINSTRUMENT: SED-Machine\nUSER: %s" % user
        header += "\nOBJECT: %s\nOUTFILE: %s" % (obj, outf)
        header += "\nOBSUTC: %s\nEXPTIME: %i" % (utc, et)
        header += "\nSKYSUB: %s" % ("On" if skysub else "Off")
        header += "\nQUALITY: %d" % qual
        header += "\nREDUCER: %s" % reducer
        header += "\nAIRMASS: %1.2f" % ec
        if 'radius_as' in ss:
            header += "\nRAPASEC: %.1f" % ss['radius_as']
        if 'xfwhm' in ss:
            header += "\nSKYLFWHMNM: %.1f" % ss['xfwhm']
        if 'yfwhm' in ss:
            header += "\nTRACEFWHMPX: %.1f" % ss['yfwhm']
        np.savetxt(outf, np.array([wl[srt], fl[srt]]).T, fmt='%8.1f  %.4e',
                   header=header)
        print("Spectrum saved to " + outf)

    if 'efficiency' in ss:
        pl.clf()
        tlab = "%s\n(Air: %1.2f | Expt: %i | Refl %d%% | Area %d cm^s)" % \
               (specname, ec, et, ss['reflectance']*100., ss['area'])
        if not no_stamp:
            pl.title(tlab)
        pl.xlabel("Wavelength [Ang]")
        pl.ylabel("SEDM efficiency (%)")
        pl.plot(ss['nm'], ss['efficiency']*100.)
        # Add drp version
        if len(drp_ver) > 0 and not no_stamp:
            ax = pl.gca()
            ax.annotate('DRP: ' + drp_ver, xy=(0.0, 0.01), xytext=(0, 0),
                        xycoords=('axes fraction', 'figure fraction'),
                        textcoords='offset points', size=6,
                        ha='center', va='bottom')
        pl.grid(True)
        pl.ioff()
        if savefig:
            outf = specname[(specname.find('_') + 1):specname.rfind('.')] + \
               '_SEDM_eff.pdf'
            pl.savefig(outf)
            print("Figure saved to " + outf)
        else:
            pl.show()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Plot extracted spectrum or data cube.

        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--cube', type=str, help='Fine correction path')
    parser.add_argument('--lambdarms', action="store_true", default=False,
                        help='Show lambda soln rms')
    parser.add_argument('--savefig', action="store_true", default=False,
                        help='Save pdf figure')
    parser.add_argument('--savespec', action="store_true", default=False,
                        help='Save spec ASCII file')
    parser.add_argument('--interact', action="store_true", default=False,
                        help='Interactively enter spectrum quality')
    parser.add_argument('--spec', type=str, help='Extracted spectrum file')
    parser.add_argument('--corrname', type=str, default='std-correction.npy')
    parser.add_argument('--redshift', type=float, default=0, help='Redshift')
    parser.add_argument('--smoothing', type=float, default=0,
                        help='Smoothing in pixels')
    parser.add_argument('--no_stamp', action="store_true", default=False,
                        help='Set to prevent plotting DRP version stamp')

    args = parser.parse_args()

    if args.cube is not None:
        check_cube(args.cube, showlamrms=args.lambdarms, savefig=args.savefig,
                   no_stamp=args.no_stamp)
    if args.spec is not None:
        check_spec(args.spec, corrname=args.corrname, redshift=args.redshift,
                   smoothing=args.smoothing,
                   savefig=args.savefig, savespec=args.savespec,
                   interact=args.interact, no_stamp=args.no_stamp)
