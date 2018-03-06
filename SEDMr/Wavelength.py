import argparse
import pdb
from multiprocessing import Pool
import numpy as np
import os
import pylab as pl
import astropy.io.fits as pf
import sys
import warnings

import datetime

import NPK.Fit
import NPK.Bar as Bar
from astropy.table import Table

from scipy.spatial import KDTree
import scipy.signal as SG

from numpy.polynomial.chebyshev import chebfit, chebval

import SEDMr.Extraction as Extraction
from scipy.interpolate import interp1d

from numpy import NaN, Inf, arange, isscalar, asarray, array


def peakdet(v, delta, x=None):
    """
    Converted from MATLAB script at http://billauer.co.il/peakdet.html
    
    Returns two arrays
    
    function [maxtab, mintab]=peakdet(v, delta, x)
    %PEAKDET Detect peaks in a vector
    %        [MAXTAB, MINTAB] = PEAKDET(V, DELTA) finds the local
    %        maxima and minima ("peaks") in the vector V.
    %        MAXTAB and MINTAB consists of two columns. Column 1
    %        contains indices in V, and column 2 the found values.
    %      
    %        With [MAXTAB, MINTAB] = PEAKDET(V, DELTA, X) the indices
    %        in MAXTAB and MINTAB are replaced with the corresponding
    %        X-values.
    %
    %        A point is considered a maximum peak if it has the maximal
    %        value, and was preceded (to the left) by a value lower by
    %        DELTA.
    
    % Eli Billauer, 3.4.05 (Explicitly not copyrighted).
    % This function is released to the public domain; Any use is allowed.
    
    """
    maxtab = []
    mintab = []

    if x is None:
        x = arange(len(v))

    v = asarray(v)

    if len(v) != len(x):
        sys.exit('Input vectors v and x must have same length')

    if not isscalar(delta):
        sys.exit('Input argument delta must be a scalar')

    if delta <= 0:
        sys.exit('Input argument delta must be positive')

    mn, mx = Inf, -Inf
    mnpos, mxpos = NaN, NaN

    lookformax = True

    for i in arange(len(v)):
        this = v[i]
        if this > mx:
            mx = this
            mxpos = x[i]
        if this < mn:
            mn = this
            mnpos = x[i]

        if lookformax:
            if this < mx-delta:
                maxtab.append((mxpos, mx))
                mn = this
                mnpos = x[i]
                lookformax = False
        else:
            if this > mn+delta:
                mintab.append((mnpos, mn))
                mx = this
                mxpos = x[i]
                lookformax = True

    return array(maxtab), array(mintab)


def read_catalog(catname):
    """Read a sextractor catalog called catname and return"""

    cat = Table.read(catname, format='ascii.sextractor')

    return cat


def hg_to_kdtree(assoc_hg_spec):
    """Convert mercury association to a KD Tree"""

    xs = []
    ys = []
    ids = []
    for id, f in enumerate(assoc_hg_spec):
        if not f.has_key(546.1): continue
        xs.append(f[546.1])
        ff = np.poly1d(f['coeff_ys'])
        ys.append(ff(xs[-1]))
        ids.append(f['seg_cnt'])

    data = np.array([xs, ys])
    return KDTree(data.T), np.array(ids)


def fiducial_wavelength():
    return 565.3


def fiducial_spectrum(lamstart=1050.0, lamratio=239./240., npx=265):
    """Return a typical SED Machine spectrum, use for interpolating on grid

    Reference:
        The equation for the fiducial spectrum looks like::

                            x
            1000 x (239/240)

    Args:
        lamstart(float): default is 1050 nm
        lamratio(float): default is 239/240, good for SEDM
        npx(int): default of 265 yields a spectrum that goes to ~ 350 nm
    """

    xx = np.arange(npx)
    cs = [1.06697198e+03,  -5.82924699e+00,   9.12740701e-03,
          -5.84149923e-06, -7.65158285e-10]
    lams = chebval(xx, cs)

    # return lamstart * lamratio**xx
    return lams


def assoc_hg_with_flats_helper(idx):
    global domedat, hgcat, guess_offset, wavetrees, n_done, update_rate

    n_done += 1
    if n_done % update_rate == 0:
        Bar.update()

    spec_pos = domedat[idx]
    if not spec_pos['ok'] and not spec_pos['bkg_ok']:
        return None
    tracefun = np.poly1d(spec_pos['coeff_ys'])
    minx = spec_pos['xs'][0]

    to_return = {}

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", category=RuntimeWarning)
        for wavelen, v in wavetrees.items():
            offset = guess_offset[wavelen]
            pt = (minx + offset, tracefun(minx+offset))
            results = v.query_ball_point(pt, 15)

            for res in results:
                x_hg, y_hg = v.data[res]
                y_trace = tracefun(x_hg)

                if np.abs(y_hg - y_trace) < 3:
                    to_return[wavelen] = x_hg

    # outstr = "\r%4.4i : %d" % (idx, len(to_return))
    # print(outstr, end="")
    # sys.stdout.flush()

    return to_return


def assoc_hg_with_flats(domedat_par, hgcat_par, guess_offset_par=None,
                        outname='assoc_Hg'):

    """Given a set of functions defining the ridgeline of dome flats
    and a list of positions of mercury lamps, provide a crude wavelength
    solution

    Args:
        guess_offset_par: Dictionary of {wavelength in nm: x position} indicating
            the rough distance from the segment start to the wavelength
    """
    global domedat, hgcat, guess_offset, wavetrees, n_done, update_rate

    domedat = domedat_par
    hgcat = hgcat_par
    if guess_offset_par is None:
        guess_offset = {365.0: 231, 404.6: 214, 435.8: 193,
                        546.1: 133, 578.00: 110}
    else:
        guess_offset = guess_offset_par

    wavetrees = {}
    for k, v in hgcat.items():
        wavetrees[k] = KDTree(v)

    print("Finding Hg lines on dome traces")
    # print("SgID   Nlines")
    n_done = 0
    update_rate = int(len(domedat) / Bar.setup(toolbar_width=74)) + 1
    p = Pool(8)
    results = p.map(assoc_hg_with_flats_helper, range(len(domedat)))
    p.close()
    Bar.done(mapped=True)
    # print("")

    for idx, spec_pos in enumerate(domedat):
        if results[idx] is None:
            continue

        for wave, pos in results[idx].items():
            spec_pos[wave] = pos

    np.save("%s.npy" % outname, [domedat])
    return spec_pos


def find_hg_spectra(lines, dYlimit=2, outname="find_spectra"):
    """Based on line ratios and positions, determine positions of Hg lines
    
    Args:
        lines: line list
        dYlimit: 

    Results:
        Writes [return value] to outname in npy format
    
    Return:
        Returns dictionary of {wavelength nm: [(pixel coords) ...], ...}
        Each wavelength has a list of pixel coords of different length.
        The coordinates are diassociated from one and other.
    """

    data = []

    for line in lines:
        data.append((line['X_IMAGE'], line['Y_IMAGE']))

    data = np.array(data)
    kdt = KDTree(data)

    reg = ""
    hg546 = []
    hg578 = []
    hg436 = []
    hg405 = []
    hg365 = []

    print("Finding Hg lines in each spectrum")
    update_rate = int(len(lines) / (Bar.setup(toolbar_width=74))) + 1
    for CNT, line in enumerate(lines):
        if CNT % update_rate == 0:
            Bar.update()

        point = [line['X_IMAGE'], line['Y_IMAGE']]
        results = kdt.query_ball_point(point, 38)
        X, Y = point
        flux = float(line['FLUX_ISO'])
        mean_flux = np.nanmedian(flux)/2.

        for resix in results:
            resX, resY = data[resix]
            resFlux = float(lines[resix]['FLUX_ISO'])
            dX = resX - X

            dY = np.abs(resY - Y)
            if dY > dYlimit:
                continue

            if resFlux < 500:
                continue

            ratio = flux/resFlux
            bright = (flux > mean_flux)
            resbright = (resFlux > mean_flux)

            if (2 < ratio < 7) and (-16 < dX < -12) and bright:
                reg += 'circle(%s,%s, 3) # color=red\n' % (X, Y)
                reg += 'circle(%s,%s, 3) # color=magenta\n' % (resX, resY)
                hg546.append((X, Y))
                hg578.append((resX, resY))
                continue

            if (1/6. < ratio < 2) and (-24 < dX < -18) and resbright:
                hg405.append((X, Y))
                hg436.append((resX, resY))
                reg += 'circle(%s,%s, 3) # color=green\n' % (X, Y)
                reg += 'circle(%s,%s, 3) # color=yellow\n' % (resX, resY)
                continue

            if (4 < ratio < 70) and (20 < dX < 38) and bright and (not resbright):
                hg365.append((resX, resY))
                reg += 'circle(%s,%s, 3) # color=blue\n' % (resX, resY)
                continue

    f = open(outname + ".reg", "w")
    f.write(reg)
    f.close()

    res = [{365.0: np.array(hg365),
            404.6: np.array(hg405),
            435.8: np.array(hg436),
            546.1: np.array(hg546),
            578.00: np.array(hg578)}]

    np.save(outname, res)
    Bar.done()

    return res[0]


def fractional_sum(FS_Y, FS_EW, FS_dat, FS_Yx1):
    """ Returns sum of FS_dat via a partial-pixel method.

    Args:
        FS_Y(float): location of trace in vertical direciton
        FS_EW(float): Width of trace around FS_Y
        FS_dat(float vector): Data to sum
        FS_Yx1(float): Start location of data

    Returns:
        Float containing the sum of data in FS_Y while taking partial
            pixels into account.

    Raises:
        Nothing.
    """

    nn = len(FS_dat)
    YB1 = FS_Y - FS_EW
    YB2 = FS_Y + FS_EW
    FSsum = 0.0
    for ii in range(nn):
        Yp1 = float(ii+FS_Yx1) - 0.5
        Yp2 = float(ii+FS_Yx1) + 0.5
        frac= 0.0
        # fully contained?
        if Yp1 > YB1 and Yp2 < YB2:
            frac = 1.0
        # pixel surrounds left boundary
        elif Yp1 < YB1 < Yp2:
            frac = (Yp2 - YB1) / (Yp2 - Yp1)
        # pixel surrounds right boundary
        elif Yp2 > YB2 > Yp1:
            frac = (YB2 - Yp1) / (Yp2 - Yp1)

        FSsum += frac * FS_dat[ii]

    return FSsum


def make_profile(sl, sigma2=4):
    """ Return a gaussian profile with the same dimensions as the slice """

    # create Gaussian profile
    profile = np.arange(np.round(sl.stop)-np.round(sl.start))
    profile = profile - (len(profile)-1)/2.0
    profile = np.exp(- profile*profile/(2*sigma2))
    # normalize it, (why mean and not sum?)
    profile /= np.mean(profile)

    return profile


def wavelength_extract_helper(SS):
    global dat, exptime, wavecalib, HDUlist, extract_width, flt_corrs, \
        n_done, update_rate

    n_done += 1
    if n_done % update_rate == 0:
        Bar.update()

    ix, flexure_x_corr_nm, flexure_y_corr_pix = SS

    ss = wavecalib[ix]
    if not ss.ok:
        return Extraction.Extraction(seg_id=ss.seg_id, ok=False,
                                     trace_sigma=ss.trace_sigma,
                                     Q_ix=ss.Q_ix, R_ix=ss.R_ix,
                                     X_as=ss.X_as, Y_as=ss.Y_as,
                                     X_pix=ss.X_pix, Y_pix=ss.Y_pix)

    minx = np.max((0, ss.xrange[0]-5))
    maxx = np.min((minx + 265, 2047))
    yfun = np.poly1d(ss.poly)

    xpos = range(minx, maxx)
    res = np.zeros(len(xpos))
    res[:] = np.nan
    resw = np.zeros(len(xpos))
    resw[:] = np.nan
    resf = np.zeros(len(xpos))
    resf[:] = np.nan
    reswf = np.zeros(len(xpos))
    reswf[:] = np.nan
    sigma2 = ss.trace_sigma * ss.trace_sigma
    if np.isnan(sigma2):
        sigma2 = 4.

    for i in range(len(xpos)):
        X = xpos[i]
        Y = yfun(X)
        if Y < 0:
            Y = 0
        if Y > 2046:
            Y = 2046
        if not np.isfinite(X) or not np.isfinite(Y):
            continue

        # Extract width requires asymmetry in the slice
        # slice[-2:3] will return elements -2 to +2 around 0
        # e.g. len(slice[-2:3]) == 5
        Ys = slice(
            int(np.max((0, np.round(Y)+flexure_y_corr_pix-extract_width))),
            int(np.min((np.round(Y)+flexure_y_corr_pix+extract_width+1, 2047))))

        # Expanded slice for fractional pixels
        Yx = slice(
            int(np.max((0, np.round(Y)+flexure_y_corr_pix-extract_width-1))),
            int(np.min((np.round(Y)+flexure_y_corr_pix+extract_width+2, 2047))))

        res[i] = np.sum(dat[Ys, X])

        profile = make_profile(Ys, sigma2=sigma2)
        resw[i] = np.sum(dat[Ys, X]*profile)
        resf[i] = fractional_sum(Y, extract_width, dat[Yx, X], Yx.start)

        profile = make_profile(Yx, sigma2=sigma2)
        reswf[i] = fractional_sum(Y, extract_width, dat[Yx, X]*profile,
                                  Yx.start)

    try:
        ll = chebval(np.arange(minx, maxx), ss.lamcoeff)
        fc = flt_corrs[ix].get_correction(ll)
    except:
        fc = 1.0
    ex = Extraction.Extraction(xrange=(minx, maxx), yrange=(yfun(xpos[0]),
                               yfun(xpos[-1])),
                               poly=ss.poly, spec=res,
                               specw=resw/fc,
                               specf=resf/fc,
                               specwf=reswf/fc,
                               seg_id=ss.seg_id, exptime=exptime, ok=True,
                               trace_sigma=ss.trace_sigma, Q_ix=ss.Q_ix,
                               R_ix=ss.R_ix, X_as=ss.X_as, Y_as=ss.Y_as,
                               X_pix=ss.X_pix, Y_pix=ss.Y_pix)

    if 'lamcoeff' in ss.__dict__ and ss.lamcoeff is not None:
        ex.lamcoeff = ss.lamcoeff.copy()
        ex.lamcoeff[0] -= flexure_x_corr_nm
        ex.lamrms = ss.lamrms
        ex.lamnrms = ss.lamnrms

    if 'mdn_coeff' in ss.__dict__ and ss.mdn_coeff is not None:
        ex.mdn_coeff = ss.mdn_coeff.copy()
        ex.mdn_coeff[0] -= flexure_x_corr_nm
        ex.lamrms = ss.lamrms
        ex.lamnrms = ss.lamnrms

    # if ex.lamrms is not None:
        # outstr = "\r%4i  %6.4f nm  %6.4f    " % (ex.seg_id, ex.lamrms, fc)
        # print(outstr, end="")
        # sys.stdout.flush()

    return ex


def wavelength_extract(HDUlist_par, wavecalib_par,
                       filename='extracted_spectra.npy',
                       flexure_x_corr_nm=0.0, flexure_y_corr_pix=0.0,
                       extract_width_par=3, inname='unknown', airmass=None,
                       flat_corrections=None):

    global dat, exptime, wavecalib, HDUlist, extract_width, flt_corrs, \
           n_done, update_rate

    extract_width = extract_width_par

    HDUlist = HDUlist_par
    wavecalib = wavecalib_par
    flt_corrs = np.copy(flat_corrections)

    dat = HDUlist[0].data
    exptime = HDUlist[0].header['EXPTIME']

    print("Applying %f nm / %f px offset" % (flexure_x_corr_nm,
                                             flexure_y_corr_pix))
    flexure_y_corr_pix = np.round(flexure_y_corr_pix)

    SSs = [(ix, flexure_x_corr_nm, flexure_y_corr_pix)
           for ix in range(len(wavecalib))]

    n_done = 0
    update_rate = int(len(SSs) / Bar.setup(toolbar_width=74)) + 1
    p = Pool(8)
    extracts = p.map(wavelength_extract_helper, SSs)
    p.close()
    Bar.done(mapped=True)

    meta = {"inname": inname, "extract_width": extract_width_par,
            "when": "%s" % datetime.datetime.now(),
            "user": os.getlogin(),
            "flexure_x_corr_nm": flexure_x_corr_nm,
            "flexure_y_corr_pix": flexure_y_corr_pix,
            "outname": filename}

    return [extracts, meta]


def extract_helper(ss):
    """ TODO: MERGE WITH wavelength_extract_helper.
    Functionallity is repeated."""
    global dat, n_done, update_rate

    n_done += 1
    if n_done % update_rate == 0:
        Bar.update()

    if not ss['ok'] and not ss['bkg_ok']:
        return Extraction.Extraction(seg_id=ss['seg_cnt'],
                                     trace_sigma=ss['trace_sigma'],
                                     bkg_ok=False, ok=False)

    minx = np.max((0, ss['xs'][0]-5))
    maxx = np.min((minx + 265, 2047))
    yfun = np.poly1d(ss['coeff_ys'])

    xpos = range(minx, maxx)
    res = np.zeros(len(xpos))
    res[:] = np.nan
    resw = np.zeros(len(xpos))
    resw[:] = np.nan
    resf = np.zeros(len(xpos))
    resf[:] = np.nan
    reswf = np.zeros(len(xpos))
    reswf[:] = np.nan
    sigma2 = ss['trace_sigma']*ss['trace_sigma']
    if np.isnan(sigma2):
        sigma2 = 4.
    extract_width = 2

    for i in range(len(xpos)):
        X = xpos[i]
        Y = yfun(X)
        if Y < 0:
            Y = 0
        if Y > 2046:
            Y = 2046
        if not np.isfinite(X) or not np.isfinite(Y):
            continue
        Ys = slice(np.max((0, np.int(Y)-3)), np.min((np.int(Y)+3, 2047)))
        # Expanded slice for fractional pixels
        Yx = slice(
            int(np.max((0, np.round(Y)-extract_width-1))),
            int(np.min((np.round(Y)+extract_width+2, 2047))))

        res[i] = np.sum(dat[Ys, X])
        # BUG: calling resw what should be resf, this is a cheap workaround
        # for now.
        resw[i] = fractional_sum(Y, extract_width, dat[Yx, X], Yx.start)

    hg_lines = {}
    for wave, pix in ss.items():
        if type(wave) is float:
            hg_lines[wave] = pix-minx-1

    return Extraction.Extraction(xrange=(minx, maxx),
                                 yrange=(yfun(xpos[0]), yfun(xpos[-1])),
                                 poly=ss['coeff_ys'], spec=res, specw=resw,
                                 hg_lines=hg_lines, seg_id=ss['seg_cnt'],
                                 bkg_ok=ss['bkg_ok'], ok=ss['ok'],
                                 trace_sigma=ss['trace_sigma'])


def extract(HDUlist, assoc_hg_spec, filename='raw_extractions'):
    global dat, n_done, update_rate

    dat = HDUlist[0].data

    n_done = 0
    update_rate = int(len(assoc_hg_spec) / Bar.setup(toolbar_width=74)) + 1
    p = Pool(8)
    extractions = p.map(extract_helper, assoc_hg_spec)
    p.close()
    Bar.done(mapped=True)

    np.save(filename, extractions)
    print("Saved extractions to %s.npy" % filename)
    return extractions


def median_fine_grid(fine, doPlot=False):
    """ Using the 7 nearest neighbors, median the wavelength solution.
    Refit the wavelength solution to this median wavelength.

    Input:
        fine: Dictionary that contains the fine wavelength solution

    Returns:
        fine.mdn_coeff contains updated chebyshev polynomial coefficients
    """

    xs = []
    ys = []
    ids = []
    for gix, g in enumerate(fine):
        # Create X, Y, and ID vector.
        if not g.ok:
            xs.append(-999)
            ys.append(-999)
            ids.append(None)
            continue

        xs.append(np.mean(g.xrange))
        ys.append(np.mean(g.yrange))
        ids.append(g.seg_id)

    xs = np.array(xs)
    ys = np.array(ys)
    ids = np.array(ids)
    dat = np.array((xs, ys)).T
    dat[dat != dat] = -999  # Correct NaNs
    KD = KDTree(dat)

    assert(len(ids) == len(fine))
    update_rate = int(len(fine) / (Bar.setup(toolbar_width=74))) + 1
    for idx, spec in enumerate(fine):
        if idx % update_rate == 0:
            Bar.update()
        seg_id = spec.seg_id
        loc = dat[seg_id-1]

        if ids[idx] is None:
            continue
        dist, nearest_ixs = KD.query(loc, k=20)

        lls = []
        num_in = 0
        for nearest in nearest_ixs[1:]:
            if nearest is None:
                continue
            if fine[nearest].hgcoef is None:
                continue
            if fine[nearest].lamcoeff is None:
                continue
            if fine[nearest].lamnrms is not None:
                if fine[nearest].lamnrms > 0.4:
                    continue
            if np.abs(fine[nearest].xrange[1] - fine[nearest].xrange[0]) < 50:
                continue

            xx = np.arange(265)
            ll = chebval(xx+fine[nearest].xrange[0], fine[nearest].lamcoeff)
            if np.any(np.diff(ll) > 0):
                continue

            lls.append(ll)
            num_in += 1
            if num_in > 30:
                break

        lls = np.array(lls, dtype=np.float32)
        if len(lls) == 0:
            spec.mdn_coeff = spec.lamcoeff
            continue

        try:
            new_lls = np.nanmedian(lls, 0)
        except:
            import pdb
            pdb.set_trace()

        diff = (lls - new_lls) / new_lls
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            stds = np.nanstd(diff, 0)
            bad = np.abs(diff/stds) > 4
        if bad.any():
            lls[bad] = np.nan
            new_lls = np.nanmedian(lls, 0)

        diff = (lls - new_lls) / new_lls
        stds = np.nanstd(diff, 0)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            bad = np.abs(diff/stds) > 2
        if bad.any():
            lls[bad] = np.nan

        new_lls = np.nanmedian(lls, 0)

        spec.mdn_coeff = chebfit(np.arange(len(new_lls))+spec.xrange[0],
                                 new_lls, 3)

    Bar.done()

    if doPlot:
        pl.figure(3)
        #plm = pl.get_current_fig_manager()
        #plm.window.wm_geometry("+670+0")
        pl.clf()
        pl.xlim(360, 550)
        pl.xlabel("Wave[nm]")
        pl.ylabel("Spec Irr [ph/10 m/nm]")
        pl.title("mdn_coeff")
        pl.ioff()
        pl.figure(4)
        pl.clf()
        pl.xlim(360, 550)
        pl.xlabel("Wave[nm]")
        pl.ylabel("Spec Irr [ph/10 m/nm]")
        pl.title("lamcoeff")
        pl.ioff()
        for g in fine:
            if g.mdn_coeff is None:
                continue
            if g.specw is None:
                continue
            if g.hgcoef is None:
                continue
            if len(g.specw) < 30:
                continue

            if g.xrange[0] < 100:
                continue
            if g.xrange[0] > 1900:
                continue
            if g.yrange[0] < 100:
                continue
            if g.yrange[0] > 1900:
                continue

            ix = np.arange(len(g.specw))
            pl.figure(3)
            ll = chebval(ix+g.xrange[0], g.mdn_coeff)
            pl.plot(ll, g.specw, '.')

            pl.figure(4)
            try:
                ll = chebval(ix+g.xrange[0], g.lamcoeff)
                pl.plot(ll, g.specw, '.')
                pl.draw()
            except:
                pass

        pl.ion()
        pl.show()

    return fine


def scale_on_547(spec):
    if spec is None:
        return None
    if spec.spec is None:
        return None
    if spec.hgcoef is None:
        return None

    ix = np.arange(0, len(spec.specw), .1)
    ix_547 = np.argmin(np.abs(chebval(ix, spec.hgcoef) - 547.0))

    scale = np.arange(len(spec.specw)) - ix[ix_547]

    return scale


def stretch_fit_helper(specno):
    """ Helper for Multiprocessing Pool"""
    global pix_coeffs, fid_ix, fs1, fs2, slopes, squares, specs, funs, n_done, update_rate

    n_done += 1
    if n_done % update_rate == 0:
        Bar.update()

    spec = specs[specno]
    if spec is None:
        return None

    xcs1 = np.zeros((len(slopes), len(squares)))
    xcs2 = np.zeros((len(slopes), len(squares)))

    s1f, s2f = funs[specno]
    for i, slope in enumerate(slopes):
        for j, square in enumerate(squares):
            # Step 4
            ns1 = s1f(fid_ix * slope + fid_ix**2 * square)
            ns2 = s2f(fid_ix * slope + fid_ix**2 * square)
            xcs1[i, j] = np.nansum(fs1*fs1*ns1*ns1)
            xcs2[i, j] = np.nansum(fs2*fs2*ns2*ns2)

    slope1, square1 = np.unravel_index(np.argmax(xcs1), xcs1.shape)
    slope2, square2 = np.unravel_index(np.argmax(xcs2), xcs2.shape)

    ix1 = fid_ix * slopes[slope1] + fid_ix**2 * squares[square1]
    ix2 = fid_ix * slopes[slope2] + fid_ix**2 * squares[square2]

    tofit = np.zeros(len(ix1))
    tofit[0:150] = ix2[0:150]
    tofit[150:] = ix1[150:]

    # Step 5
    fc = chebfit(fid_ix, tofit, 2)
    xc = (xcs1[slope1, square1] + xcs2[slope2, square2]) / np.nansum(fs1)

    # outstr = "\r%4.0i = %1.3e : %1.3f %+1.2e, %1.3f %+1.2e" % (specno, xc,slopes[slope1], squares[square1], slopes[slope2], squares[square2])
    # print(outstr,)
    # sys.stdout.flush()

    if False:
        pl.figure(4)
        pl.plot(fid_ix, s1f(tofit)+s2f(tofit))

        import pdb
        pdb.set_trace()

    return fc, xc


def get_stretched(fid_ix, coeffs, spec1, spec2=None):
    """ Re-interpolate spec1 and spec2 onto the fiducial index set by
    stretching and scaling the spectra

    This is used after a call that looks something like:
        ll = np.load("Hg_ext.npy")
        lx = np.load("Xe_ext.npy")
        gg = rough_grid(ll)
        coeffs, fix = stretch_set(gg, lx)

        s1, s2 = get_stretched(fix, coeffs[5], ll[5], lx[5])
    
    Returns:
        s1, s2 [float(len(fid_ix))]: Spectra
    """

    newix = chebval(fid_ix, coeffs)
    ix = scale_on_547(spec1)
    sf1 = interp1d(ix, spec1.specw, bounds_error=False, fill_value=0)
    if spec2 is None:
        sf2 = lambda x: 0.0
    else:
        sf2 = interp1d(ix, spec2.specw, bounds_error=False, fill_value=0)

    s1 = sf1(newix)
    s2 = sf2(newix)

    return s1/s1.mean()*spec1.specw.mean(), s2/s2.mean()*spec2.specw.mean()


def stretch_set(Hg_set, Xe_set, mscales=None):
    """ Shift and stretch spectra so they have the same pixel index 

    Steps:
    1. Shift spectra without interpolation onto a common pixel grid with
     index 0 being a prominent Hg line.
    2. Create a set of interpolating functions for each of the spectra
     in the set. These functions are called as flux(pixel).
    3. Brute force search for the best stretch and 2nd-order polynomial
     coeff for each spectrum in the set, compared to the fiducial spectrum.
    4. Brute force will measure the above for the xenon and mercury spectra
     indepdentnly.
    5. Then we stitch together a master pixel shift based on a chebyshev fit
     of the xenon and mercury spectra. The chebyshev fit creates a function
     such that f(pixel) --> pixel. On this new grid, all spectra are close
     to the same.
    """

    # Global is for interprocess comms
    global pix_coeffs, fid_ix, fs1, fs2, slopes, squares, specs, funs, n_done, update_rate

    assert(len(Hg_set) == len(Xe_set))

    # Step 1
    gridded_set = []
    for index, hg in enumerate(Hg_set):
        gridded_set.append(None)
        xe = Xe_set[index]
        if hg.spec is None:
            continue
        if xe.spec is None:
            continue
        if hg.hgcoef is None:
            continue

        assert(len(hg.spec) == len(xe.spec))

        ix = np.arange(len(hg.specw))
        ix_547 = np.argmin(np.abs(chebval(ix, hg.hgcoef) - 547.0))
        if ix_547 > 180:
            continue
        if ix_547 < 100:
            continue

        scale = np.arange(len(hg.spec)) - ix_547

        gridded_set[-1] = (scale, hg.spec, xe.spec)

    fid_index = int(len(gridded_set)/2)
    fiducial = gridded_set[fid_index]
    while fid_index < len(gridded_set) and (fiducial is None or
                                                    len(fiducial[0]) != 265):
        fid_index += 1
        fiducial = gridded_set[fid_index]

    if fiducial is None:
        print("No good fiducial")
        quit()
    else:
        print("Fiducial index: %d" % fid_index)

    print("Number of extractions: %d" % len(gridded_set))

    # Step 2
    pl.figure(3)
    pl.clf()

    def interp_functions():
        funs = []
        for num, element in enumerate(gridded_set):
            funs.append(None)
            if element is None:
                continue

            ix, s1, s2 = element
            s1f = interp1d(ix, s1, fill_value=np.nan, bounds_error=False)
            s2f = interp1d(ix, s2, fill_value=np.nan, bounds_error=False)

            funs[-1] = (s1f, s2f)

        return funs

    funs = interp_functions()

    pix_coeffs = []
    fid_ix, fs1, fs2 = fiducial
    # Note 0.001 over 250 pixels is about a quarter pixel
    # max of 0.5e-3 over 150**2 is a handful of pixels.
    slopes = np.arange(0.95, 1/.95, .001)
    squares = np.linspace(-1e-3, 1e-3, 15)

    specs = gridded_set
    update_rate = int(len(specs) / Bar.setup(toolbar_width=74)) + 1
    p = Pool(8)
    results = p.map(stretch_fit_helper, range(len(specs)))
    p.close()
    Bar.done(mapped=True)

    pix_coeffs = []
    xcs = []
    for res in results:
        if res is None:
            pix_coeffs.append(None)
            xcs.append(None)
        else:
            pix_coeffs.append(res[0])
            xcs.append(res[1])

    return pix_coeffs, fiducial, xcs


# TODO -- REFIT Wavelengths given a good guess

# Linelist for calib lamps (currently not using He)
linelist = {
    "He": [587.5, 667.8, 706.5],
    "Cd": [467.8, 479.9, 508.5],
    "Hg": [578, 546.1, 435.8, 404.6, 365],
    "Xe": [764, 828]
}


def snap_solution_into_place(PARS):
    """ Return new Chebyshev coefficients for best fit wavelength solution """

    global n_done, update_rate

    n_done += 1
    if n_done % update_rate == 0:
        Bar.update()

    if PARS is None:
        return None, None, None
    ixs, coef, lamp_spec = PARS

    # Use 5-parameter Gaussian to fit line peaks
    fitfun = NPK.Fit.gaussian5
    resfun = NPK.Fit.mpfit_residuals(fitfun)

    def fitlinelist(cc, printout=False):
        """ Fit a line list using input coefficients as a first guess at
        dispersion function. 
        
        INPUTS:
        
        cc - Chebyshev coeffs to use
        
        KEYWORDS:
        
        printout - set to True to get running status on line-fitting
        """

        # Initial solution
        lams = chebval(ixs, cc)

        # array for results
        results = []

        # Loop over lamps: Hg, Xe, Cd, He?
        for lampname, lampspec in lamp_spec.items():

            # Loop over lines for the given lamp
            for line in linelist[lampname]:

                # Find the appropriate line
                lix = np.argmin(np.abs(line-lams))

                # Avoid ends
                if lix < 5 or lix > (len(ixs)-5):
                    continue

                # Get sub-array of spectrum
                cutout = slice(lix-4, lix+4)
                xct = ixs[cutout]               # Position (x)
                sct = lampspec[cutout]          # Flux (y)

                # Initial guess for fit
                parguess = [
                    {'value': sct.max()-sct.min()},     # Scale
                    {'value': xct[sct.argmax()]},       # Centroid!!
                    {'value': 1.2},                     # Sigma
                    {'value': sct.min()},               # offset
                    {'value': 0}]                       # slope

                # Do the fitting
                pars = NPK.Fit.mpfit_do(resfun, xct, sct, parguess,
                                        error=np.sqrt(np.abs(sct)))
                # Skip bad fits
                if pars.status != 1:
                    continue

                # Gather results for each line in the spectrum
                results.append((line, pars.params[1], pars.perror[1]))

                # Print results
                # if printout:
                #   out_str = \
                #       "\rLine: %7.2f nm, Xpos: %8.3f px, dXpos: %6.3f px" % \
                #         (line, pars.params[1], pars.perror[1])
                #   print(out_str, end="")
                #   sys.stdout.flush()
            # END: for line in linelist[lampname]:
        # END: for lampname, lampspec in lamp_spec.items():

        # Convert results into an numpy array
        results = np.array(results)

        # Skip if we only solved for fewer than 3 lines
        if len(results) < 3:
            return cc, np.nan, np.nan

        # Break out the parameters
        LS = results[:, 0]       # Reference wavelength
        IXS = results[:, 1]      # Pixel position of line centroid
        WS = results[:, 2]       # Error in position in pixels

        # Calculate number of coefficients to use in fit
        nc = len(LS)-1
        if nc > 5:
            nc = 4

        # Fit dispersion function
        newcoef = chebfit(IXS, LS, nc, w=WS)

        # Residuals
        res = np.abs(chebval(IXS, newcoef) - LS)
        nres = res/LS

        # Error checking on 365 nm line (not currently used)
        # if nres[LS==365] < 200:
        #    sl = np.where(LS==365)[0]
        #    np.delete(IXS, sl)
        #    np.delete(LS, sl)
        #    np.delete(WS, sl)

        #    newcoef = chebfit(IXS, LS, nc, w=WS)
        #    res = np.abs(chebval(IXS, newcoef) - LS)
        #    nres = res/LS

        return newcoef, nres, res
    # END: def fitlinelist(cc, printout=False):

    # Iterate three times on fitting dispersion function
    newcoef, nres, res = fitlinelist(coef)
    newcoef, newnres, newres = fitlinelist(newcoef)
    newcoef, newnres, newres = fitlinelist(newcoef, printout=True)

    # return the coefs and the RMS (normalized and unnormalized)
    return newcoef, np.sqrt(np.mean(newnres*newnres)), \
           np.sqrt(np.mean(newres*newres))


def snap_solution_into_place_all(fine, Hgs, Xes, Cds=None, Hes=None):
    """ Get a global solution for dispersion function
    using Hg, Xe, possibly Cd data
    """

    global n_done, update_rate

    # Data structure for passing to fitting function
    PARS = []

    # Loop over extractions in fine
    for i in range(len(fine)):

        # Add an empty record
        PARS.append(None)

        # Skip bad extractions
        if fine[i].xrange is None:
            continue
        if fine[i].mdn_coeff is None:
            continue

        # Get X values
        ixs = np.arange(*fine[i].xrange)
        # Get starting coeffs
        coef = fine[i].mdn_coeff
        # Start with Hg and Xe spectra
        lamp_spec = {"Hg": fine[i].specw, "Xe": Xes[i].specw}

        # Add Cd and/or He if we have them
        if Cds is not None:
            lamp_spec["Cd"] = Cds[i].specw
        if Hes is not None:
            lamp_spec["He"] = Hes[i].specw

        # Populate with X-values, coefficients, and spectra
        PARS[-1] = (ixs, coef, lamp_spec)
    # END: for i in range(len(fine)):

    n_done = 0
    update_rate = int(len(PARS) / Bar.setup(toolbar_width=74)) + 1
    # Map onto thread pool
    p = Pool(8)
    results = p.map(snap_solution_into_place, PARS)
    p.close()
    Bar.done(mapped=True)
    # print("")

    # Calculate diagnostics of fits
    min_rms = 1.e9
    max_rms = 0.
    good_rms = []

    # loop over fits
    for i, res in enumerate(results):

        # Re-populate fine cube
        fine[i].lamcoeff = res[0]       # Fit coefficients
        fine[i].lamnrms = res[1]        # Normalized RMS
        fine[i].lamrms = res[2]         # Un-normalized RMS (nm)

        # If we have a good fit
        if fine[i].lamrms is not None:

            # Compile residuals
            good_rms.append(res[2])
            if res[2] > max_rms:
                max_rms = res[2]
            if min_rms > res[2] > 0.:
                min_rms = res[2]
    # END: for i, res in enumerate(results):

    # Calculated and report RMS statistics
    avg_rms = np.nanmean(good_rms)
    print("<RMS>: %8.3f nm, RMS(min): %8.3f nm, RMS(max): %8.3f nm" %
          (float(avg_rms), min_rms, max_rms))
    print("")

    return fine


def fit_he_lines(SS, guesses=None, Ncutout=7):
    if guesses is None:
        return fit_known_lines(SS, {587.5: -18, 667.8: -48, 706.5: -60},
                               Ncutout)
    else:
        return fit_known_lines(SS, guesses, Ncutout)


def fit_cd_lines(SS, guesses=None, Ncutout=5):
    if guesses is None:
        return fit_known_lines(SS, {467.8: 41, 479.9: 33, 508.5: 18}, Ncutout)
    else:
        return fit_known_lines(SS, guesses, Ncutout)


def fit_hg_lines(SS, guesses=None, Ncutout=7):
    if guesses is None:
        return fit_known_lines(SS, {578: -14, 546.1: 0, 435.8: 60, 404.6: 81},
                               Ncutout)
    else:
        return fit_known_lines(SS, guesses, Ncutout)


def fit_xe_lines(SS, guesses=None, Ncutout=5):
    if guesses is None:
        return fit_known_lines(SS, {764: -77, 828: -93}, Ncutout)
    else:
        return fit_known_lines(SS, guesses, Ncutout)


def fit_known_lines(SS, guesses, Ncutout):
    """ Fit line positions based on guess pixel positions against a fiducial spectrum 

    This function is mapable
    
    Args:
        SS is a list of two elements containting:
            spix(int[Ns]): Index values (pixel) of spec
            spec(float[Ns]): Flux from corresponding index values
        guesses: Dictionary containing {wavelength in nm: pixel guess position}
    Returns:
        [ (wavelength in nm, centroid position in pixel) ] with length 
            len(guesses.keys()).

        Intent is for this list to be fit with polynomials to construct the
        wavelength solution.
        
    """

    if SS is None:
        return None

    # Break out parameters
    spix, spec = SS

    # Array for results
    res = []

    # Use 5-parameter Gaussian to fit line peaks
    fitfun = NPK.Fit.gaussian5
    resfun = NPK.Fit.mpfit_residuals(fitfun)

    # Loop over lines to fit
    for lam, guesspos in guesses.items():

        # Get sub-array of line to fit
        cutout = np.where(np.abs(spix - guesspos) < Ncutout)[0]
        xct = spix[cutout]      # Line position (X)
        sct = spec[cutout]      # Line flux (Y)

        # Initial guess at fit parameters
        parguess = [
            {'value': sct.max()-sct.min()},     # Scale
            {'value': xct[sct.argmax()]},       # Centroid!!
            {'value': 1.3},                     # Sigma
            {'value': sct.min()},               # offset
            {'value': 0}]                       # slope

        # Do the peak fitting
        pars = NPK.Fit.mpfit_do(resfun, xct, sct, parguess)

        # Skip bad fits
        if pars.status != 1:
            continue

        # Store results: ref wavelength, and centroid
        res.append([lam, pars.params[1]])
    # END: for lam, guesspos in guesses.items():

    # Convert results into numpy array
    res = np.array(res)

    return res


def fit_all_lines(fiducial, hg_spec, xe_spec, xxs, cd_spec=None, he_spec=None):
    """Fit mercury + xenon lines to a set of spectra that are put on a common fiducial grid.

    Args:
        fiducial(int[Nf]): Fiducial index values range from about -130 to + 130
        hg_spec(Extraction[Ne]): Mercury spectra
        xe_spec(Extraction[Ne]): Xenon spectra
        xxs(Gridded[Ne]): Results from stretch_set() call

    """

    global n_done, update_rate

    # Verify that all extractions share the same dimensionality
    assert(len(hg_spec) > 500)
    assert(len(hg_spec) == len(xe_spec))
    if cd_spec is not None:
        assert(len(hg_spec) == len(cd_spec))
    if he_spec is not None:
        assert(len(hg_spec) == len(he_spec))

    assert(len(xxs) == len(hg_spec))

    # Arrays in which to accumulate extracted spectra
    Hgs = []
    Xes = []
    Cds = []
    Hes = []

    print("Stretching spectra")

    # Loop over extracted spectra
    for i in range(len(hg_spec)):

        # Placeholders for bad extractions
        if hg_spec[i] is None or xxs[i] is None:
            Hgs.append(None)
            Xes.append(None)
            Cds.append(None)
            Hes.append(None)
            continue

        # Stretch using Hg and Xe
        s1, s2 = get_stretched(fiducial, xxs[i], hg_spec[i],
                               spec2=xe_spec[i])

        # Stretch using Hg and Cd
        if cd_spec is not None:
            s1, s3 = get_stretched(fiducial, xxs[i], hg_spec[i],
                                   spec2=cd_spec[i])
            # Store stretched Cd spectra
            Cds.append((fiducial, s3))

        # Stretch using Hg and He
        if he_spec is not None:
            s1, s4 = get_stretched(fiducial, xxs[i], hg_spec[i],
                                   spec2=he_spec[i])
            # Store stretched He spectra
            Hes.append((fiducial, s4))

        # Store stretched Hg, Xe spectra
        Hgs.append((fiducial, s1))
        Xes.append((fiducial, s2))

    # Now find all the line centroids in stretched spectra

    # Start with Hg
    print("Getting Hg line positions")
    n_done = 0
    update_rate = int(len(Hgs) / Bar.setup(toolbar_width=74)) + 1
    p = Pool(8)
    hg_locs = p.map(fit_hg_lines, Hgs)
    p.close()
    Bar.done(mapped=True)

    # Now get Xe lines
    print("Getting Xe line positions")
    n_done = 0
    update_rate = int(len(Xes) / Bar.setup(toolbar_width=74)) + 1
    p = Pool(8)
    xe_locs = p.map(fit_xe_lines, Xes)
    p.close()
    Bar.done(mapped=True)

    # Cd next, if present
    if cd_spec is not None:
        print("Getting Cd line positions")
        n_done = 0
        update_rate = int(len(Cds) / Bar.setup(toolbar_width=74)) + 1
        p = Pool(8)
        cd_locs = p.map(fit_cd_lines, Cds)
        p.close()
        Bar.done(mapped=True)

    # Finally He, if present
    if he_spec is not None:
        print("Getting He line positions")
        n_done = 0
        update_rate = int(len(Hes) / Bar.setup(toolbar_width=74)) + 1
        p = Pool(8)
        he_locs = p.map(fit_he_lines, Hes)
        p.close()
        Bar.done(mapped=True)

    # Now fit all lines for dispersion function
    print("Determining best fit to all lines")

    # Store results
    fits = []
    residuals = []
    rss = []

    # Accumulate fit statistics
    mean_rms = 0.
    max_rms = 0.
    min_rms = 1.e9

    # Loop over spectra
    for i in range(len(hg_spec)):

        # Placeholders for results
        fits.append(None)
        rss.append(None)

        # Get Hg, Xe positions
        hg = hg_locs[i]
        xe = xe_locs[i]

        # Ancillary positions, if supplied
        if cd_spec is not None:
            cd = cd_locs[i]
        if he_spec is not None:
            he = he_locs[i]

        # Check to be sure our data are good
        if hg is None:
            continue
        if xe is None:
            continue
        if cd_spec is not None and cd is None:
            continue
        if he_spec is not None and he is None:
            continue

        # Create fitting input vectors
        try:

            # Start with Hg, Xe
            ll = np.concatenate([hg[:, 0], xe[:, 0]])     # Wavelengths
            ix = np.concatenate([hg[:, 1], xe[:, 1]])     # Positions

            # Add Cd and He if supplied
            if cd_spec is not None:
                ll = np.concatenate([ll, cd[:, 0]])
                ix = np.concatenate([ix, cd[:, 1]])
            if he_spec is not None:
                ll = np.concatenate([ll, he[:, 0]])
                ix = np.concatenate([ix, he[:, 1]])

        # Skip if we can't concatenate
        except:
            continue

        # Get fitting weights
        weights = np.ones(len(ll))
        weights[ll < 400] = 0.03          # Lower weight for extreme blue
        weights[ll > 860] = 0.1           # Lower weight for extreme red

        # Calculate number of coefficients for fit
        nc = 3
        if len(ll) < 5:
            nc = 2

        # Fit dispersion
        ff = chebfit(ix, ll, nc, w=weights)

        # Store fit coefficients
        fits[-1] = ff

        # Calculate residuals
        res = {}
        for jx, l in enumerate(ll):
            res[l] = chebval(ix[jx], ff) - l
        residuals.append(res)

        rss[-1] = np.sqrt(np.sum((chebval(ix, ff) - ll)**2))
        mean_rms += rss[-1]
        if rss[-1] < min_rms:
            min_rms = rss[-1]
        if rss[-1] > max_rms:
            max_rms = rss[-1]

        if False:
            pl.figure(3)
            pl.step(chebval(fiducial, ff), Hgs[i][1]+Hes[i][1]+Xes[i][1])
            pl.figure(3)
            pl.step(chebval(fiducial, ff), Hgs[i][1]+Hes[i][1]+Xes[i][1])

    pl.clf()
    XS = []
    YS = []
    CS = []
    for i in range(len(hg_spec)):
        try:
            x = np.mean(hg_spec[i].xrange)
            y = np.mean(hg_spec[i].yrange)
            c = rss[i]
        except:
            XS.append(None)
            YS.append(None)
            continue

        XS.append(x)
        YS.append(y)
        CS.append(c)

    mean_rms /= 1. * len(hg_spec)
    print("<RMS>: %8.3f nm, RMS(min): %8.3f nm, RMS(max): %8.3f nm" %
          (mean_rms, min_rms, max_rms))
    print("")

    return np.array(fits), residuals, np.array(rss), (XS, YS)


def coeffs_to_spec(fix, gridded, rgrd_coef, lam_coef):
    """ Returns the spectrum given the unadulterated spectrum and fits.

    Args:
        fix(array) -- Fiducial index positions
        gridded(Extraction) -- Extracted spectrum to plot
        rgrd_coef(list) -- Chebyshev polynomial coefficients that represent
            the stretch of the spectrum to put onto the fiducial index
            array (fix)
        lam_coef(list) -- Chebyshev polynomial coefficients that convert
            fix to wavelength in nm

    Returns:
        Tuple containing the wavelength and spectrum
        (wavelength, spectrum) this is in the units of the fit coefficients.
        
    """

    gix = scale_on_547(gridded)
    gsp = gridded.specw
    newix = chebval(fix, rgrd_coef)

    CC = chebfit(newix, fix, 3)
    XX = chebval(gix, CC)
    LL = chebval(XX, lam_coef)

    return LL, gsp


def plot_grid(grid, Xes=None):

    pl.figure(2)
    pl.clf()
    pl.figure(1)
    pl.clf()
    pl.ioff()
    ixs = []
    for index, g in enumerate(grid):
        if index < 800:
            continue
        if index > 900:
            continue
        if g.specw is None:
            continue
        if len(g.specw) == 0:
            continue
        if g.hgcoef is None:
            continue
        ix = np.arange(len(g.specw))
        ix_547 = np.argmin(np.abs(chebval(ix, g.hgcoef) - 547))

        if ix_547 > 180:
            continue
        if ix_547 < 100:
            continue

        ixs.append(ix_547)

        ix = np.arange(len(g.specw)) - ix_547

        add = 0
        if Xes is not None:
            try:
                add = Xes[index].specw
            except:
                pass
            if add is None:
                add = 0
        pl.figure(2)
        pl.plot(ix, g.specw+add, 'x')

    pl.show()

    pl.figure(3)
    pl.clf()
    pl.ion()
    pl.hist(ixs, 50)


def rough_grid_helper(ix):
    global extractions, lines, n_ext, tot_std, min_std, max_std

    ex = extractions[ix]

    if not ex.ok and not ex.bkg_ok:
        return

    xs = []
    ys = []
    for line in lines:
        if line in ex.hg_lines:
            X = ex.hg_lines[line]
            if X + ex.xrange[0] > 2048:
                ex.ok = False
                continue
            xs.append(ex.hg_lines[line])
            ys.append(line)

    xs = np.array(xs)
    ys = np.array(ys)

    if len(xs) == 0:
        # print(ex.xrange, ex.yrange, ex.hg_lines)
        return

    with warnings.catch_warnings():
        warnings.simplefilter("ignore", np.polynomial.polyutils.RankWarning)
        coef = chebfit(xs, ys, 2)
        vals = chebval(xs, coef)

    ex.hgcoef = coef
    err = (vals - ys)

    rough_std = np.std(err)
    tot_std += rough_std
    n_ext += 1
    if rough_std > max_std:
        max_std = rough_std

    if rough_std < min_std:
        min_std = rough_std

    # outstr = "\r%4.4i: rough RMS: %6.3f nm" % (ex.seg_id, np.std(err))
    # print(outstr,)
    # sys.stdout.flush()


def rough_grid(data, the_lines=None, outname='rough_wavelength.npy'):
    """Shift the extractions onto a coarse common pixel grid"""

    global extractions, lines, n_ext, tot_std, min_std, max_std

    extractions = data
    if the_lines is None:
        lines = [365.0, 404.6, 435.8, 546.1, 578]
    else:
        lines = the_lines
    tot_std = 0.
    min_std = 1.e9
    max_std = 0.
    n_ext = 0

    for i in range(len(extractions)):
        rough_grid_helper(i)

    if n_ext > 0:
        avg_std = tot_std / (1. * n_ext)
    else:
        print("Warning - no good extractions")
        avg_std = 0.
    print("<STD>: %6.3f nm, STD(min): %6.3f nm, STD(max): %6.3f nm" %
          (avg_std, min_std, max_std))
    print("")

    try:
        np.save(outname, extractions)
        print("Wrote %s.npy" % outname)
    except:
        pass
    return extractions


def RMS(vec):
    return np.sqrt(np.sum((vec-np.mean(vec))**2))


def fit_spectra_Hg_Xe(Hgs, Xes, kdtree, kdseg_ids, plot=False,
                      outname='fit_spectra'):

    assert(len(Hgs) == len(Xes))

    if plot:
        pl.figure(1)
        pl.clf()
        pl.figure(2)
        pl.clf()
        pl.figure(1)

    ixs = np.arange(0, 260, .01)
    for spec_ix in range(len(Hgs)):
        hg = Hgs[spec_ix]
        xe = Xes[spec_ix]
        # seg_id = hg.seg_id

        coeff = hg.hgcoef
        if coeff is None:
            print(spec_ix)
            continue

        prev_rms = np.inf
        rms = np.inf
        for i in range(25):

            offsets = measure_offsets(hg.specw+xe.specw, coeff, plot=False)

            # keys = np.array(offsets.keys())
            # vals = np.array(offsets.values())
            # print(np.sort(keys+vals))

            lams = chebval(ixs, coeff)
            pixs = []
            meas = []

            for lam, off in offsets.items():
                if (lam+off) > 1000: continue
                if (lam+off) < 340: continue
                ix = ixs[np.argmin(np.abs(lam-lams))]
                pixs.append(ix)
                meas.append(lam+off)

            if np.any(np.isfinite(meas) == False):
                break
            if np.any(pixs == 0):
                break
            if len(pixs) < 2:
                break

            if i < 2:
                coeff = chebfit(pixs, meas, 3)
            else:
                coeff = chebfit(pixs, meas, 4)
            diff = chebval(pixs, coeff) - meas
            rms = RMS(diff)
            # print("%4.0i %6.3f" % (spec_ix, rms))

            if not np.isfinite(rms):
                pdb.set_trace()

            if rms < 0.35:
                break
            if np.abs(rms - prev_rms) < 0.02:
                break
            prev_rms = rms

        if plot:
            pl.figure(2)
            lam = chebval(np.arange(len(hg.specw)), coeff)
            spec = hg.specw+xe.specw
            pl.plot(lam, spec)
            pl.figure(1)

        Hgs[spec_ix].lamcoeff = coeff
        Hgs[spec_ix].lamrms = rms
        print("-- %4.0i %6.3f" % (spec_ix, rms))

    if plot:
        pl.show()

    np.save(outname, Hgs)
    print("Wrote %s.npy" % outname)
    return Hgs


def measure_offsets(spec, coeff, plot=False):
    """Measure wavelength offsets of Hg and Xe extractions

    First uses a crude peak finding code to identify the locations of the 
    Xe lines.

    Then, the two Xe complexes are fit individually
        
    Returns:
        A dictionary of {line_wavelength_nm: offset_nm}.
    """

    preffun = lambda x: ((x[1])**2 + (x[0]-120)**2)/1e4 + 1.
    resfun830 = NPK.Fit.mpfit_residuals(xe_830nm, preffun=preffun)

    # Preference function is designed to bias against adjusting the 
    preffun = lambda x: (np.abs(x[0]-100)/100 + np.abs(x[1])/30 +
                         np.abs(x[4]-100)/100)/3e4 + 1.0
    resfun890 = NPK.Fit.mpfit_residuals(xe_890nm, preffun=preffun)
    resfung = NPK.Fit.mpfit_residuals(NPK.Fit.gaussian5)

    pk_pxs = SG.argrelmax(spec, order=10)[0]
    pk_lam = chebval(pk_pxs, coeff)
    ix = np.arange(len(spec))
    ll = chebval(ix, coeff)

    if plot:
        pl.figure(3)
        pl.clf()
        pl.plot(ll, spec)
        for p in pk_lam:
            print("-", p)
            pl.axvline(p,color='blue')

    offsets = {}

    # resfung parameters:

    for index, lam in enumerate([365.0, 404.6, 435.8, 546.1, 578]):
        roi = (ll > (lam-9)) & (ll < (lam+9))

        toll = ll[roi]
        tofit = spec[roi]
        if len(toll) == 0:
            continue
        offset = np.min(tofit)
        slope = 0
        sigma = 2.5
        scale = np.max(tofit)-offset

        parinfo = [{'value': scale, 'limited': [1, 0], 'limits':[0, 0]},
                   {'value': lam, 'limited': [1, 1], 'limits':[lam-15, lam+15]},
                   {'value': sigma, 'limited': [1, 1], 'limits': [1.5, 3.5]},
                   {'value': offset, 'limited': [1, 0], 'limits': [0, 0]},
                   {'value': slope}]
        res = NPK.Fit.mpfit_do(resfung, toll, tofit, parinfo)

        if res.status > 0:
            offsets[lam] = lam-res.params[1]

    guess_lam = pk_lam[1]
    ok = ((guess_lam-35) < ll) & (ll < (guess_lam+45))

    if np.any(ok):
        sigma = 80
        lamoffset = guess_lam-830
        offset = np.min(spec[ok])
        peak = np.max(spec[ok]) - offset

        parinfo = [{'value': sigma, 'limited': [1, 1], 'limits': [25, 250]},
                   {'value': lamoffset, 'limited': [1, 1], 'limits': [-50, 50]},
                   {'value': offset, 'limited': [1, 0], 'limits': [0, 0]},
                   {'value': peak, 'limited': [1, 0], 'limits': [0, 0]}]

        res = NPK.Fit.mpfit_do(resfun830, ll[ok], spec[ok], parinfo)
        thell = ll[ok]
        ps = res.params.copy()
        ps[1] = 0
        thespec = xe_830nm(ps, thell)
        cl = np.sum(thespec*thell)/np.sum(thespec)
        offsets[cl] = -res.params[1]

        if plot:
            pl.plot(ll[ok], xe_830nm(res.params, ll[ok]), 'x-')

    guess_lam = pk_lam[0]
    ok = ((guess_lam-55) < ll) & (ll < (guess_lam+100))

    if np.any(ok):
        sigma = 100
        lamoffset = guess_lam-890
        if lamoffset > 5:
            lamoffset = 5
        if lamoffset < -5:
            lamoffset = -5
        offset = np.min(spec[ok])
        peak = np.max(spec[ok]) - offset
        p930 = 500

        parinfo = [{'value': sigma, 'limited': [1, 1], 'limits': [25, 260]},
                   {'value': lamoffset, 'limited': [1, 1], 'limits': [-50, 50]},
                   {'value': offset, 'limited': [1, 0], 'limits': [0, 0]},
                   {'value': peak, 'limited': [1, 0], 'limits': [0, 0]},
                   {'value': p930, 'limited': [1, 0], 'limits': [0, 0]}]

        res = NPK.Fit.mpfit_do(resfun890, ll[ok], spec[ok], parinfo)

        thell = ll[ok]
        ps = res.params.copy()
        thespec = xe_890nm(ps, thell)
        cl = np.sum(thespec*thell)/np.sum(thespec)
        offsets[cl] = -res.params[1]

        if plot:
            print(res.status)
            print(res.niter, res.perror)
            print(res.params)
            print(res.debug)
            print(res.errmsg)
            pl.plot(ll[ok], xe_890nm(res.params, ll[ok]), 'x-')

    if plot:
        pdb.set_trace()
    return offsets


def read_spec_loc(spec_loc_fname):
    """Returns the structure for the location of each spectrum"""

    return np.load(spec_loc_fname)


def xe_830nm(p, lam):
    """Xe complex near 830 nm.

    See: http://www2.keck.hawaii.edu/inst/lris/arc_calibrations.html
    """

    sigma, lamoffset, offset, peak = p

    sig = (2000*np.exp(-(lam-820.63-lamoffset)**2/sigma) +
           1000*np.exp(-(lam-828.01-lamoffset)**2/sigma) +
           100*np.exp(-(lam-834.68-lamoffset)**2/sigma) +
           180*np.exp(-(lam-840.82-lamoffset)**2/sigma))
    sig = sig/max(sig)*peak + offset

    return sig


def save_fitted_ds9(fitted, outname='fine'):

    ds9 = 'physical\n'
    for ix,S in enumerate(fitted):
        if S.lamcoeff is None: continue
        if S.xrange is None: continue

        xvals = np.arange(*S.xrange)

        res = chebval(xvals, S.lamcoeff)
        invcoeffs = chebfit(res, xvals, 4)
        pxs = chebval([400., 500., 600., 656.3, 700., 800., 900.], invcoeffs)

        for ix, px in enumerate(pxs):
            X = px
            Y = np.poly1d(S.poly)(X)

            if ix == 3:
                ds9 += 'point(%s,%s) # point=cross text={%s} color=blue\n' % \
                    (X, Y, S.seg_id)
            elif ix == 0:
                ds9 += 'point(%s,%s) # point=box color=red\n' % (X, Y)
            elif ix == 6:
                ds9 += 'point(%s,%s) # point=circle color=red\n' % (X, Y)
            else:
                ds9 += 'point(%s,%s) # point=cross color=red\n' % (X ,Y)

    f = open(outname+".reg", "w")
    f.write(ds9)
    f.close()


def xe_890nm(p, lam):
    """Xe complex near 890 nm.

    See: http://www2.keck.hawaii.edu/inst/lris/arc_calibrations.html
    """

    sigma, lamoffset, offset, peak, p937 = p

    sig = (5000*np.exp(-(lam-881.94-lamoffset)**2/sigma) +
           1000*np.exp(-(lam-895.22-lamoffset)**2/sigma) +
           1000*np.exp(-(lam-904.54-lamoffset)**2/sigma) +
           1900*np.exp(-(lam-916.26-lamoffset)**2/sigma) +
           p937*np.exp(-(lam-937.42-lamoffset)**2/sigma))
    sig = sig/max(sig)*peak + offset

    return sig


def assign_fit_to_spectra(target, gridded, rss, fix, stretch, lamfit):
    """ Put wavelength solution into target
    
    Args:
        target(list of Extraciton): Target list of extractions
        gridded(list of Extractions[Ne]): Gridded spectra
        rss (list of float[Ne]): RSS wavelength fit error
        fix (list of float): Fiducial spectral index  for "stretch"
        stretch (list of float [Ne x Nfloat]): Float polynomials to stretch
        spectra onto
        lamfit (list of float [Ne x Ncoeff]): Wavelength fitting function
        coefficients
        
    Returns:
        New grid with coefficients assigned to it

    """

    for i in range(len(gridded)):

        try:
            L, S = coeffs_to_spec(fix, gridded[i], stretch[i], lamfit[i])
        except:
            continue

        cc = chebfit(np.arange(*target[i].xrange), L, 5)
        target[i].lamcoeff = cc
        target[i].lamrms = rss[i]

    return target


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
    """Wavelength.py performs:

        1. Rough wavelength calibration based on Mercury lamp lines.
        2. Fine wavelength calbiration based on Xenon and Mercury lamp lines.
        3. Extraction of spectra with a fine wavelength calibration.

        For each step, various report files are written, often ds9 region 
        files.

        the _coarse.npy file contains a length-1 list with a dictionary. The
            dictionary contains {wavelength: [(XY point)]} with all XY points
            for a given Hg emission line wavelength.

        the assoc_Hg.npy file contains a length-1 list with the associated
            Hg solutions from the above *_coarse file. 

        --dome comes from FindSpectra.py

        for rough set --dome [npy], --hgcat [txt], and --outname 
        for fine set --xefits [fits], --hefits [fits] --cdfits [fits] --hgassoc [npy], and --outname
    """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('step', type=str, help='One of [rough|fine|extract]')
    parser.add_argument('--infile', type=str,
                        help='Infile name, purpose depends on step')
    parser.add_argument('--dome', type=str,
                        help='Dome segment definition npy file. Used in rough.')
    parser.add_argument('--hgfits', type=str, help='Name of mercury fits file')
    parser.add_argument('--xefits', type=str, help='Name of xenon fits file')
    parser.add_argument('--cdfits', type=str, help='Name of cadmium fits file')
    parser.add_argument('--hefits', type=str, help='Name of helium fits file')
    parser.add_argument('--hgcat', type=str,
                        help='Name of mercury sextractor catalog')
    parser.add_argument('--coarse', type=str,
                        help='Name of coarse Hg solution [npy]')
    parser.add_argument('--hgassoc', type=str,
                        help='Name of coarse Hg solution [npy]')
    parser.add_argument('--fine', type=str, help='Name of fine solution [npy]')
    parser.add_argument('--toextract', type=str,
                        help='Name of fine solution [npy]')
    parser.add_argument('--outname', type=str,
                        help='Prefix name of output file')

    args = parser.parse_args()
    outname = args.outname

    if args.step == 'rough':
        print("\nROUGH WAVELENGTH SOLUTION\n")
        hgfits = args.hgfits
        catname = args.hgcat
        catalog = read_catalog(catname)
        spec_loc_fname = args.dome
        spec_loc = read_spec_loc(spec_loc_fname)
        hg_spec = find_hg_spectra(catalog, outname=outname)
        assoc_hg_with_flats(spec_loc, hg_spec)
    elif args.step == 'fine':
        print("\nFINE WAVELENGTH SOLUTION\n")
        XeDat = pf.open(args.xefits)
        HgDat = pf.open(args.hgfits)
        if args.cdfits is not None:
            CdDat = pf.open(args.cdfits)
        else:
            CdDat = None
        if args.hefits is not None:
            HeDat = pf.open(args.hefits)
        else:
            HeDat = None
        assoc_hg_spec = np.load(args.hgassoc)[0]

        print("Extracting Hg")
        Hg_E = extract(HgDat, assoc_hg_spec, filename="Hg_ext_"+args.outname)

        print("Extracting Xe")
        Xe_E = extract(XeDat, assoc_hg_spec, filename="Xe_ext_"+args.outname)

        Cd_E = None
        if CdDat is not None:
            print("Extracting Cd")
            Cd_E = extract(CdDat, assoc_hg_spec, filename="Cd_ext_"+args.outname)
        He_E = None
        if HeDat is not None:
            print("Extracting He")
            He_E = extract(HeDat, assoc_hg_spec, filename="He_ext_"+args.outname)

        print("Getting rough grid")
        gridded = rough_grid(Hg_E)

        print("Getting stretched set for Xe")
        stretchset, fiducial, xcors = stretch_set(gridded, Xe_E)
        fix, f1, f2 = fiducial
        print("Fitting all lines")
        fits, residuals, rss, locs = fit_all_lines(fix,
                                                   gridded,
                                                   Xe_E,
                                                   stretchset,
                                                   cd_spec=Cd_E,
                                                   he_spec=He_E)

        print("Assigning fit to spectra")
        result = assign_fit_to_spectra(Hg_E, gridded, rss,
                                       fix, stretchset, fits)

        print("Getting median fine grid")
        result = median_fine_grid(result)

        print("Snapping solution into place")
        snap_solution_into_place_all(result, Hg_E, Xe_E, Cds=Cd_E, Hes=He_E)
        np.save(args.outname, result)
        print("Wrote %s.npy" % args.outname)

    elif args.step == 'extract':

        print("\nEXTRACTION\n")
        fitted = np.load(args.fine)
        hdu = pf.open(args.toextract)
        outname = args.outname

        ww = wavelength_extract(hdu, fitted, filename=outname,
            inname=args.toextract)
        ww[1]['airmass'] = hdu.header['airmass']
        np.save(outname, ww)
        print("Wrote %s.npy" % outname)

    sys.exit()

