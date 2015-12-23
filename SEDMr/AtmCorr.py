
import argparse, os, sys
import numpy as np
import scipy.stats
import scipy.signal
import pylab as pl
import pyfits as pf
import datetime
import os
import sets
import warnings

from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.interpolate import interp1d

import Wavelength
import NPK.Standards as SS

def handle_create(outname=None, filelist=[], plot_filt=False):
    ''' Create standard star correction. Units are erg/s/cm2/Ang '''

    if outname is None: outname='atm-corr.npy'

    ll = Wavelength.fiducial_spectrum()
    corrs = []
    corr_vals =[]
    legend=["corr",]
    filt_legend=["orig"]
    maxnm = 915.
    for file in filelist:
        ''' Filename is sp_STD-nnnnn_obs*.npy '''

        # Try to read input file
        try: data = np.load(file)[0]
        except:
            raise Exception("Not passed a spectrum file, the file should start with sp_")

        # Check for calculated correction
        if "std-correction" not in data.keys():
            print "No std-correction vector in %s" % file
            continue
        correction = data['std-correction']

        # What was the maximum wavelength?
        if "std-maxnm" in data.keys():
            if data['std-maxnm'] > maxnm: maxnm = data['std-maxnm']

        # Convert file named sp_STD-nnnn_obs* to a name compaitable
        # with Standards.py
        pred = file.lstrip("sp_STD-")
        pred = pred.split("_")[0]
        legend.append(pred)
        pred = pred.lower().replace("+","").replace("-","_")
        print file, pred, pred in SS.Standards

        if pred not in SS.Standards:
            print "File named '%s' is reduced to '%s' and no such standard seems to exist."  % (file, pred)
            continue

        # Record median correction in ROI
        ROI = (ll > 600) & (ll < 850)
        corr_vals.append(np.median(correction[ROI]))
        
        # Normalize each correction by the median value
        correction /= corr_vals[-1]
        corrs.append(correction)
        
    # END: for file in filelist:

    # Rescale each correction and scale to erg/s/cm2/Ang
    corrs = np.array(corrs)
    erg_s_cm2_ang = corrs * np.median(corr_vals) * 1e-16

    # Take the median of the correction vectors
    the_corr = scipy.stats.nanmedian(erg_s_cm2_ang,0)
    
    # Fit red end unless we are calibrated out there
    if not np.isfinite(the_corr).all() and maxnm < 1000.0:
        print "Fitting red end"
        # Fit polynomial to extrapolate correction
        redend = (ll > 880) & (ll < maxnm)
        ss = np.poly1d(np.polyfit(ll[redend], the_corr[redend], 2))

        redend = (ll>maxnm)
        the_corr[redend] = ss(ll[redend])

    # Plot data, if requested
    if plot_filt:
        pl.figure(2)
        pl.clf()
        pl.grid(True)
        pl.ylim(1e-20,1e-15)
        pl.semilogy(ll, the_corr, linewidth=1)

        pl.xlabel("Wavelength [nm]")
        pl.ylabel("Correction [erg/s/cm cm/Ang]")
        pl.title("Correct ph/10 m/nm to erg/s/cm2/Ang")

    # Now clean over the Balmer series
    balmers = [656.3, 486.1, 434.0, 410.2, 397.0]
    #balmers = [656.3, 486.1, 434.0]
    eps = 0.02
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
        the_corr[to_fix] = fit(ll[to_fix])
        #pl.plot(ll[to_fix], the_corr[to_fix])
        #pl.show()

    if plot_filt:
        pl.semilogy(ll, the_corr*4., linewidth=2)
        filt_legend.append("Balmer*4")

    # Filter correction to remove spectral features and leave response alone
    the_corr = scipy.signal.savgol_filter(the_corr, 9, 5)

    if plot_filt:
        pl.semilogy(ll, the_corr*2., linewidth=2)
        filt_legend.append("Filtered*2")
        pl.semilogy(ll, the_corr, linewidth=2)
        filt_legend.append("Filtered")
        pl.legend(filt_legend)
        pl.show()

    # Plot data
    pl.figure(1)
    pl.clf()
    pl.grid(True)
    pl.ylim(1e-20,1e-15)
    pl.semilogy(ll, the_corr, linewidth=4)
    for ix,e in enumerate(erg_s_cm2_ang):
        pl.semilogy(ll, e*corr_vals[ix]/np.mean(corr_vals))

    pl.xlabel("Wavelength [nm]")
    pl.ylabel("Correction [erg/s/cm cm/Ang]")
    pl.title("Correct ph/10 m/nm to erg/s/cm2/Ang")
    pl.legend(legend)

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        pl.savefig("Standard_Correction.pdf")

    print "Mean cor: %10.3g, Sigma cor: %10.3g" % (np.mean(corr_vals) * 1e-16, np.std(corr_vals)*1e-16)
    maxnm = np.min([maxnm,np.max(ll)])
    print "Max nm: %7.2f" % maxnm

    # Construct result
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        res = {"nm": ll,
            "maxnm": maxnm,
            "correction": the_corr,
            "doc": "Correct ph/10 m/nm to erg/2/cm2/ang",
            "Nspec": len(corrs),
            "correction_std": np.nanstd(erg_s_cm2_ang,0),
            "outname": outname,
            "files": filelist,
            "when": '%s' % datetime.datetime.now(),
            "user": os.getlogin()
            }

    np.save(outname, [res])
    return res

def handle_summary(outname=None, filelist=[]):
    if outname is None: outname = 'correction.npy'

    keepers = []
    for file in filelist:
        print file
        f = np.load(file)[0]
        if np.nanmin(f['std-correction']) < 9:
            keepers.append(f)
            print "Keeping %s" % file


    corl = keepers[0]['nm'].copy()
    cor = np.zeros((len(corl), len(keepers)))
    cor[:,0] = keepers[0]['std-correction'].copy()
    for ix, keeper in enumerate(keepers[1:]):
        f = interp1d(keeper['nm'], keeper['std-correction'], bounds_error=False,
            fill_value = np.nan)

        cor[:,ix] = f(corl)

    cs = np.nanmean(cor,1)
    ccs = chebfit(corl, cs, 6)

    cor = [{"nm": corl, "cor": cs, "coeff": ccs}]
    np.save(outname, cor)

def handle_corr(filename, outname='std-correction.npy', objname=None) :

    if outname is None: outname = "corr_" + filename
    dat = np.load("sp_" + filename)[0]
    erg_s_cm2_ang = dat['std-correction'] * 1e-16
    maxnm = dat['std-maxnm']

    object = filename.split('-')[1].split('.')[0]

    pl.figure(1)
    pl.clf()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        pl.xlim(np.nanmin(dat['nm'])-10., maxnm + 10.)
        pl.semilogy(dat['nm'], erg_s_cm2_ang)
        pl.ylim(1e-20,1e-15)
        pl.grid(True)

        pl.xlabel("Wavelength [nm]")
        pl.ylabel("Correction [erg/s/cm cm/Ang]")
        pl.title("ph/10 m/nm to erg/s/cm2/Ang from %s" % object)

        pl.savefig("corr_" + filename.rstrip(".npy") + ".pdf")

    # Construct result
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        res = {"nm": dat['nm'],
            "maxnm": maxnm,
            "correction": erg_s_cm2_ang,
            "doc": "Correct ph/10 m/nm to erg/2/cm2/ang",
            "Nspec": 1,
            "correction_std": np.nanstd(erg_s_cm2_ang),
            "outname": outname,
            "files": filename,
            "when": '%s' % datetime.datetime.now(),
            "user": os.getlogin()
            }

    np.save(outname, [res])


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=\
        '''


        ''', formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('process', type=str, help='Process [CORR|SUM|CREATE]')
    parser.add_argument('--A', type=str, help='FITS A file')
    parser.add_argument('--outname', type=str, help='Prefix output name')
    parser.add_argument('--std', type=str, help='Name of standard')
    parser.add_argument('--files', type=str, metavar='file', nargs='+',
        help='Name of standard')

    args = parser.parse_args()


    if args.process == 'CORR':
        # Take atmospheric correction out and store in a separate file
        handle_corr(args.A, outname=args.outname, objname=args.std)

    if args.process == 'SUM':
        handle_summary(outname=args.outname, filelist=args.files)

    if args.process == 'CREATE':
        # Create the atmospheric correction.
        handle_create(outname=args.outname, filelist=args.files)

