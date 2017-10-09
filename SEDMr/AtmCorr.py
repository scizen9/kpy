"""Process standard star corrections

 The three options for processing are:
    * CORR    - extract a single correction vector
    * SUM     - obsolete version of CREATE
    * CREATE  - create an ensemble correction

Functions:
    * :func:`handle_corr` output a single correction
    * :func:`handle_create` create a calibration spectrum
    * :func:`handle_summary` obsolete version of handle_create

Note:
    This is used as a python script as follows::

        usage: AtmCorr.py [-h] [--A A] [--outname OUTNAME] [--std STD]
                  [--files file [file ...]]
                  process

        positional arguments:
          process               Process [CORR|SUM|CREATE]

        optional arguments:
          -h, --help            show this help message and exit
          --A A                 FITS file to correct
          --outname OUTNAME     Prefix output name
          --std STD             Name of standard
          --files file [file ...]
                                list of spectrum files: sp_STD-*.npy

"""
import argparse
import datetime
import os
import warnings

import numpy as np
import pylab as pl
import scipy.signal
from numpy.polynomial.chebyshev import chebfit
from scipy.interpolate import interp1d

import NPK.Standards as Stds


def handle_corr(filename, outname='corrected.npy'):
    """Output single correction. """

    if outname is None:
        outname = "corr_" + filename
    dat = np.load("sp_" + filename)[0]
    if "std-correction" not in dat.keys():
        print("Not a known standard extraction, returning")
        return
    erg_s_cm2_ang = dat['std-correction']
    maxnm = dat['std-maxnm']

    objname = filename.split('-')[1].split('.')[0]

    pl.figure(1)
    pl.clf()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        pl.xlim(np.nanmin(dat['nm']) - 10., maxnm + 10.)
        pl.semilogy(dat['nm'], erg_s_cm2_ang)
        pl.ylim(1e-20, 1e-15)
        pl.grid(True)

        pl.xlabel("Wavelength [nm]")
        pl.ylabel("Correction [erg/s/cm cm/Ang]")
        pl.title("ph/10 m/nm to erg/s/cm2/Ang from %s" % objname)

        pl.savefig("corr_" + os.path.splitext(filename)[0] + ".pdf")

    # Construct result
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        res = {
            "nm": dat['nm'],
            "maxnm": maxnm,
            "correction": erg_s_cm2_ang,
            "doc": "Correct ph/10 m/nm to erg/s/cm2/ang",
            "Nspec": 1,
            "correction_std": np.nanstd(erg_s_cm2_ang),
            "outname": outname,
            "files": filename,
            "when": '%s' % datetime.datetime.now(),
            "user": os.getlogin()
            }

    np.save(outname, [res])
    print("Wrote %s.npy" % outname)


def handle_create(outname=None, filelist=None, plot_filt=False):
    """Create standard star correction. Units are erg/s/cm2/Ang

    Read in the std-correction vector for each standard star, massage
    the correction (see below) and output an ensemble correction.

    Args:
        outname (str): name for resulting correction array (def: atm-corr.npy)
        filelist (list): list of standard star extractions (e.g. sp_STD-\*.npy)
        plot_filt (bool): set to plot intermediate filtered steps (def: False)

    Returns:
        dict: A dictionary containing the ensemble correction and some
        details::

            {   "nm": wavelengths (nm),
                "maxnm": maximum wavelength,
                "correction": the correction spectrum,
                "doc": "Correct ph/10 m/nm to erg/s/cm2/ang",
                "Nspec": number of spectra used,
                "correction_std": STDev of ensemble corrections,
                "outname": outname (see Parameters),
                "files": filelist (see Parameters),
                "when": timestamp string,
                "user": user who ran the process    }

    Note:
        Writes a \*.npy file with the resulting correction

    """

    if outname is None:
        outname = 'atm-corr.npy'

    ll = None
    corrs = []
    corr_vals = []
    legend = ["corr", ]
    filt_legend = ["orig"]
    maxnm = 915.
    drp_ver = None
    if filelist is None:
        filelist = []
    for ifile in filelist:
        """ Filename is sp_STD-nnnnn_obs*.npy """

        # Try to read input file
        try:
            data = np.load(ifile)[0]
        except:
            print("Not able to load %s. " % ifile)
            continue

        # Check quality of extraction
        if data["quality"] > 0:
            print("Bad std extraction in %s" % ifile)
            continue

        # Check for calculated correction
        if "std-correction" not in data.keys():
            print("No std-correction vector in %s" % ifile)
            continue

        if ll is None:
            ll = data['nm']
            correction = data['std-correction']
        else:
            corf = interp1d(data['nm'], data['std-correction'],
                            bounds_error=False, fill_value=np.nan)
            correction = corf(ll)

        # What was the maximum wavelength?
        if "std-maxnm" in data.keys():
            if data['std-maxnm'] > maxnm:
                maxnm = data['std-maxnm']
        # Convert file named sp_STD-nnnn_obs* to a name compatable
        # with Standards.py.  First strip of "sp_STD-"
        pred = ifile[7:]
        pred = pred.split("_")[0]
        legend.append(pred)
        pred = pred.lower().replace("+", "").replace("-", "_")
        print(ifile, pred, pred in Stds.Standards)
        # Are we in our list of standard stars?
        if pred not in Stds.Standards:
            print("File named '%s' is reduced to '%s' and no such standard "
                  "seems to exist." % (ifile, pred))
            continue
        # get DRP version
        if drp_ver is None:
            drp_ver = data['drp_version']
        # Record median correction in ROI
        roi = (ll > 600) & (ll < 850)
        corr_vals.append(np.nanmedian(correction[roi]))
        # Normalize each correction by the median value
        correction /= corr_vals[-1]
        corrs.append(correction)
    # END: for ifile in filelist:

    if len(corrs) > 0:
        print("Fitting %d standards" % len(corrs))

        # Rescale each correction and scale to erg/s/cm2/Ang
        corrs = np.array(corrs)
        erg_s_cm2_ang = corrs * np.nanmedian(corr_vals)
        # Take the median of the correction vectors
        the_corr = np.nanmedian(erg_s_cm2_ang, 0)
        # Fit red end unless we are calibrated out there
        if not np.isfinite(the_corr).all() and maxnm < 1000.0:
            print("Fitting red end")
            # Fit polynomial to extrapolate correction
            redend = (ll > 880) & (ll < maxnm) & np.isfinite(the_corr)
            ss = np.poly1d(np.polyfit(ll[redend], the_corr[redend], 2))
            # Insert extrapolation back into correction vector
            redend = (ll > maxnm)
            the_corr[redend] = ss(ll[redend])

        # Plot data, if requested
        if plot_filt:
            pl.figure(2)
            pl.clf()
            pl.grid(True)
            pl.ylim(1e-20, 1e-15)
            pl.semilogy(ll, the_corr, linewidth=1)
            pl.xlabel("Wavelength [nm]")
            pl.ylabel("Correction [erg/s/cm cm/Ang]")
            pl.title("Correct ph/10 m/nm to erg/s/cm2/Ang")

        """
        # Now clean over the Balmer series (skip this for now)
        balmers = [656.3, 486.1, 434.0, 410.2, 397.0]
        #balmers = [656.3, 486.1, 434.0]
        for balmer in balmers:
            #pl.figure(3)
            line_ROI = sets.Set(np.where(np.abs((ll-balmer)/balmer) < 0.01)[0])
            broad_ROI = sets.Set(np.where(np.abs((ll-balmer)/balmer) < 0.04)[0])
            broad_ll = list(broad_ROI)
            broad_ll.sort()
            #pl.plot(ll[broad_ll], the_corr[broad_ll])
            around_line_ROI = list(broad_ROI - line_ROI)
            #pl.plot(ll[around_line_ROI], the_corr[around_line_ROI],'^')
            fit = np.poly1d(np.polyfit(ll[around_line_ROI],
                the_corr[around_line_ROI], 5))
            to_fix = list(line_ROI)
            #the_corr[to_fix] = fit(ll[to_fix])
            #pl.plot(ll[to_fix], the_corr[to_fix])
            #pl.show()
            """

        # Plot intermediate correction
        if plot_filt:
            pl.semilogy(ll, the_corr * 4., linewidth=2)
            filt_legend.append("Balmer*4")
        # Filter nans
        fin = np.isfinite(the_corr)
        # Filter correction to remove spectral features and leave response alone
        the_corr[fin] = scipy.signal.savgol_filter(the_corr[fin], 9, 5)
        # Plot intermediate correction
        if plot_filt:
            pl.semilogy(ll, the_corr * 2., linewidth=2)
            filt_legend.append("Filtered*2")
            pl.semilogy(ll, the_corr, linewidth=2)
            filt_legend.append("Filtered")
            pl.legend(filt_legend)
            pl.show()

        # Plot data
        pl.figure(1)
        pl.clf()
        pl.grid(True)
        pl.ylim(1e-20, 1e-15)
        pl.semilogy(ll, the_corr, linewidth=4)
        for ix, e in enumerate(erg_s_cm2_ang):
            pl.semilogy(ll, e * corr_vals[ix] / np.mean(corr_vals))

        pl.xlabel("Wavelength [nm]")
        pl.ylabel("Correction [erg/s/cm cm/Ang]")
        pl.title("Correct ph/10 m/nm to erg/s/cm2/Ang")
        if drp_ver is not None:
            ax = pl.gca()
            ax.annotate('DRP: ' + drp_ver, xy=(0.0, 0.01), xytext=(0, 0),
                        xycoords=('axes fraction', 'figure fraction'),
                        textcoords='offset points', size=6,
                        ha='center', va='bottom')
        pl.legend(legend)

        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            pl.savefig("Standard_Correction.pdf")

        print("Mean cor: %10.3g, Sigma cor: %10.3g" %
              (float(np.mean(corr_vals)),
               float(np.std(corr_vals))))
        maxnm = np.min([maxnm, np.max(ll)])
        print("Max nm: %7.2f" % maxnm)

        # Construct result
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            res = {
                "nm": ll,
                "maxnm": maxnm,
                "correction": the_corr,
                "doc": "Correct ph/10 m/nm to erg/s/cm2/ang",
                "Nspec": len(corrs),
                "correction_std": np.nanstd(erg_s_cm2_ang, 0),
                "outname": outname,
                "files": filelist,
                "when": '%s' % datetime.datetime.now(),
                "user": os.getlogin()
            }

        np.save(outname, [res])
        print("Wrote %s.npy" % outname)
        return res
    else:
        print("No good sp_STD*.npy files, nothing generated")


def handle_summary(outname=None, filelist=None):
    """Extract all std-correction vectors and create ensemble (obsolete)"""

    if outname is None:
        outname = 'correction.npy'

    if filelist is not None:
        # Get a list of good correction vectors
        keepers = []
        for ifile in filelist:
            print(ifile)
            f = np.load(ifile)[0]
            # Correction values should be small
            # if "std-correction" not in data.keys():
            if np.nanmin(f['std-correction']) < 9:
                keepers.append(f)
                print("Keeping %s" % ifile)

        # Set fiducial wavelengths from first correction vector
        corl = keepers[0]['nm'].copy()
        # Insert first vector
        cor = np.zeros((len(corl), len(keepers)))
        cor[:, 0] = keepers[0]['std-correction'].copy()
        # Insert the rest of the vectors
        for ix, keeper in enumerate(keepers[1:]):
            f = interp1d(keeper['nm'], keeper['std-correction'],
                         bounds_error=False, fill_value=np.nan)
            cor[:, ix] = f(corl)
        # Create mean correction
        cs = np.nanmean(cor, 1)
        # Fit wl coefficients
        ccs = chebfit(corl, cs, 6)
        # Create output
        cor = [{"nm": corl, "cor": cs, "coeff": ccs}]
        np.save(outname, cor)
        print("Wrote %s.npy" % outname)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Process standard star corrections in one of three ways:
  CORR    - extract a single correction vector
  SUM     - obsolete version of CREATE
  CREATE  - create an ensemble correction
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('process', type=str, help='Process [CORR|SUM|CREATE]')
    parser.add_argument('--A', type=str, help='FITS file to correct')
    parser.add_argument('--outname', type=str, help='Prefix output name')
    parser.add_argument('--std', type=str, help='Name of standard')
    parser.add_argument('--files', type=str, metavar='file', nargs='+',
                        help='list of spectrum files: sp_STD-*.npy')

    args = parser.parse_args()

    if args.process == 'CORR':
        # Take atmospheric correction out and store in a separate file
        handle_corr(args.A, outname=args.outname)

    if args.process == 'SUM':
        handle_summary(outname=args.outname, filelist=args.files)

    if args.process == 'CREATE':
        # Create the atmospheric correction.
        handle_create(outname=args.outname, filelist=args.files)
