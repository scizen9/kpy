
import argparse
import os
import glob
import numpy as np
import pylab as pl
import pyfits as pf
import itertools
import warnings

from astropy.coordinates import Angle
from astropy.time import Time
from numpy.polynomial.chebyshev import chebfit, chebval
from scipy.interpolate import interp1d
from scipy.ndimage import filters
import SEDMr.Wavelength as Wavelength
import SEDMr.Spectra as SedSpec
import SEDMr.GUI as GUI
import NPK.Util
import NPK.Standards as Stds
import NPK.Atmosphere as Atm

# Nadia imports
from scipy.interpolate import griddata
import scipy.optimize as opt


def reject_outliers(data, m=2.):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    d[d != d] = 100000.
    s = d/mdev if mdev else 0.
    return data[s < m]


def get_ellipse_xys(ell):
    a = ell[0]
    b = ell[1]

    pts = np.zeros((361, 2))
    beta = -ell[4] * np.pi / 180.
    sin_beta = np.sin(beta)
    cos_beta = np.cos(beta)
    alpha = np.radians(np.r_[0.:360.:1j * 361])

    sin_alpha = np.sin(alpha)
    cos_alpha = np.cos(alpha)

    pts[:, 0] = ell[2] + (a * cos_alpha * cos_beta - b * sin_alpha * sin_beta)
    pts[:, 1] = ell[3] + (a * cos_alpha * sin_beta + b * sin_alpha * cos_beta)

    return pts


def atm_dispersion_positions(prlltc, pos, leff, airmass):
    """ Return list of (X,Y) positions indicating trace of
    atmospheric dispersion

    Args:
        prlltc (float): parralactic angle in Angle class
        pos (float,float): (x,y) position of source in arcsec at wave leff
        leff (float): Effective wavelength, micron
        airmass (float): Atmospheric airmass, airmass=1 means no dispersion

    Returns:
        (list): List of source positions in arcsec: [ (x0,y0) ... (xn,yn) ]

        Note: if airmass=1, then the list is equivalent of [ pos ]

    """
    print "LamEff %.1f microns, Airmass %.3f" % (leff, airmass)

    # Compute the AD for two pieces of the spectrum:
    # from 0.38 microns to leff and from leff to 0.95 microns
    blue_ad = NPK.Util.atm_disper(0.38, leff, airmass)
    red_ad = NPK.Util.atm_disper(0.95, leff, airmass)
    print 'Blue AD is %1.1f", Red AD is %1.1f" PRLLTC %3.1f' % (blue_ad, red_ad,
                                                                prlltc.degree)

    dx = -np.sin(prlltc.radian)
    dy = np.cos(prlltc.radian)

    delta = 0.1
    bpos = np.array(pos) - np.array([dx, dy]) * blue_ad  # * delta

    positions = []
    nstep = np.int(np.round((blue_ad - red_ad)/delta))
    if nstep == 0:
        nstep = 1
    for step in xrange(nstep):
        t = [bpos[0] + step * dx * delta, bpos[1] + step * dy * delta]
        positions.append(t)

    dx = positions[0][0] - positions[-1][0]
    dy = positions[0][1] - positions[-1][1]

    print "DX %2.1f, DY %2.1f, D %2.1f" % (dx, dy, np.sqrt(dx*dx + dy*dy))
    return positions


def gaussian_2d(xdata_tuple, amplitude, xo, yo,
                sigma_x, sigma_y, theta, offset):
    """
    Produces a 2D gaussian centered in xo, yo with the parameters specified.
    xdata_tuple: coordinates of the points where the 2D Gaussian is computed.

    """
    (x, y) = xdata_tuple
    xo = float(xo)
    yo = float(yo)
    a = (np.cos(theta)**2)/(2*sigma_x**2) + (np.sin(theta)**2)/(2*sigma_y**2)
    b = -(np.sin(2*theta))/(4*sigma_x**2) + (np.sin(2*theta))/(4*sigma_y**2)
    c = (np.sin(theta)**2)/(2*sigma_x**2) + (np.cos(theta)**2)/(2*sigma_y**2)
    g = offset + amplitude*np.exp(-(a * ((x-xo)**2) +
                                    2 * b * (x-xo) * (y-yo) +
                                    c * ((y-yo)**2)))
    return g.ravel()


def identify_spectra_gauss_fit(spectra, prlltc=None, lmin=400., lmax=900.,
                               airmass=1.0, sigfac=3.0):
    """ 
    Returns index of spectra picked by Guassian fit.
    
    NOTE: Index is counted against the array, not seg_id
    """

    status = 0

    pl.ioff()

    kt = SedSpec.Spectra(spectra)

    # Get X,Y positions (arcsec) and summed values between lmin and lmax
    xs, ys, vs = kt.to_xyv(lmin=lmin, lmax=lmax)

    xi = np.linspace(np.nanmin(xs), np.nanmax(xs), 100)
    yi = np.linspace(np.nanmin(ys), np.nanmax(ys), 200)

    x, y = np.mgrid[np.nanmin(xs):np.nanmax(xs):100j,
                    np.nanmin(ys):np.nanmax(ys):200j]

    points = zip(xs, ys)
    values = vs

    grid_vs = griddata(points, values, (x, y), method='linear')
    grid_vs[np.isnan(grid_vs)] = np.nanmean(grid_vs)
    print("grid_vs min, max, mean: %f, %f, %f" %
          (np.nanmin(grid_vs), np.nanmax(grid_vs), np.nanmean(grid_vs)))

    # pl.plot(np.nansum(grid_vs, axis=1))
    # pl.plot(np.nansum(grid_vs, axis=0))
    # pl.show()
    # Initialize the first guess for the Gaussian
    xo = xi[np.argmax(np.nansum(grid_vs, axis=1))]
    yo = yi[np.argmax(np.nansum(grid_vs, axis=0))]
    sigma_x = 1.3
    sigma_y = 1.
    amplitude = np.nanmax(grid_vs)
    print("initial guess: z,a,b,x,y,theta: %f, %f, %f, %f, %f, %f" %
          (amplitude, sigma_x, sigma_y, xo, yo, 0.))

    # create data
    initial_guess = (amplitude, xo, yo, sigma_x, sigma_y, 0,
                     np.nanmean(grid_vs))

    try:
        popt, pcov = opt.curve_fit(gaussian_2d, (x, y),
                                   grid_vs.flatten(), p0=initial_guess)
    except RuntimeError:
        print "ERROR: unable to fit Gaussian"
        print "Using initial guess"
        status = 3
        popt = initial_guess

    xc = popt[1]
    yc = popt[2]
    if xc < -30. or xc > 30. or yc < -30. or yc > 30.:
        print "ERROR: X,Y out of bounds: %f, %f" % (xc, yc)
        print "Using initial position: %f, %f" % (xo, yo)
        xc = xo
        yc = yo
        status = 1
    pos = (xc, yc)

    # get 3-sigma extent
    a = popt[3]*sigfac
    b = popt[4]*sigfac
    if a > 30. or b > 30.:
        print "ERROR: A,B out of bounds: %f, %f" % (a, b)
        print "Using initial values: %f, %f" % (sigma_x, sigma_y)
        a = sigma_x * sigfac
        b = sigma_y * sigfac
        status = 2
    theta = popt[5]
    z = popt[0]

    # report position and shape
    ellipse = (a, b, xc, yc, theta * (180. / np.pi))
    print("PSF FIT on IFU:  z,a,b,x,y,theta = %f, %f, %f, %f, %f, %f" %
          (z, a, b, xc, yc, theta*180./np.pi))

    leffmic = (lmax+lmin)/2000.0    # convert to microns

    if prlltc is not None:
        positions = atm_dispersion_positions(prlltc, pos, leffmic, airmass)
    else:
        positions = [pos]

    all_kix = []
    for the_pos in positions:
        all_kix.append(list(find_positions_ellipse(kt.KT.data, xc, yc, a, b,
                                                   -theta)))

    all_kix = list(itertools.chain(*all_kix))
    kix = list(set(all_kix))
    print "found this many spaxels: %d" % len(kix)

    return kt.good_positions[kix], pos, positions, ellipse, status


def find_positions_ellipse(xy, h, k, a, b, theta):
    """
    xy: Vector with pairs [[x0, y0], [x1, y1]] of coordinates.
    a: semi-major axis of ellipse in X axis.
    b: semi-minor axis of ellipse in Y axis.
    h: central point ellipse in X axis.
    k: central point ellipse in Y axis.
    theta: angle of rotation of ellipse in radians (clockwise).
    """
    positions = np.arange(len(xy))
    x = xy[:, 0]
    y = xy[:, 1]
    dist = ((x-h) * np.cos(theta) + (y-k) * np.sin(theta)) ** 2 / (a ** 2) + \
           ((x-h) * np.sin(theta) - (y-k) * np.cos(theta)) ** 2 / (b ** 2)

    return positions[dist < 1]


def identify_spectra_gui(spectra, radius=2., scaled=False, bgd_sub=True,
                         lmin=650., lmax=700., cmin=-300, cmax=300, prlltc=None,
                         objname=None, airmass=1.0, nosky=False, message=None):
    """ Returns index of spectra picked in GUI.

    NOTE: Index is counted against the array, not seg_id
    """

    # Set ellipse parameters
    elfl = glob.glob("ell_STD-*.npy")
    if len(elfl) > 0:
        inell = np.load(elfl[0])
        inell[1] = radius*(inell[1]/inell[0])
        inell[0] = radius
        inell[2] = 0.
        inell[3] = 0.
        print "Loaded ellipse parameters from %s" % elfl[0]
    else:
        inell = (radius, radius*(3./5.), 0., 0., 20.5)
        print "Using default ellipse parameters"
    print "\nStarting with a %s arcsec semimajor axis" % radius

    # Get spectral extractions
    kt = SedSpec.Spectra(spectra)

    # Scale cube
    if not scaled:
        s = GUI.ScaleCube(kt, bgd_sub=bgd_sub, lmin=lmin, lmax=lmax,
                          objname=objname)
        scaled = s.scaled

        if scaled:
            cmin = s.cmin
            cmax = s.cmax

    # Print message if needed
    if message is not None:
        print message

    # Get positions
    g = GUI.PositionPicker(kt, bgd_sub=bgd_sub, ellipse=inell, scaled=scaled,
                           lmin=lmin, lmax=lmax, cmin=cmin, cmax=cmax,
                           objname=objname, nosky=nosky)
    pos = g.picked
    nosky = g.nosky
    ellipse = g.ellipse

    print "Final semimajor axis (arcsec) = %4.1f" % ellipse[0]
    if nosky:
        print "Sky subtraction off"
    else:
        print "Sky subtraction on"

    leffmic = (lmax+lmin)/2000.0    # Convert to microns

    if prlltc is not None:
        positions = atm_dispersion_positions(prlltc, pos, leffmic, airmass)
    else:
        positions = [pos]

    all_kix = []
    for the_pos in positions:
        all_kix.append(list(find_positions_ellipse(kt.KT.data,
                                                   ellipse[2], ellipse[3],
                                                   ellipse[0], ellipse[1],
                                                   -ellipse[4]*(np.pi/180.))))

    all_kix = list(itertools.chain(*all_kix))
    kix = list(set(all_kix))

    stats = {"nosky": nosky, "scaled": scaled,
             "lmin": lmin, "lmax": lmax,
             "cmin": cmin, "cmax": cmax}

    return kt.good_positions[kix], pos, positions, ellipse, stats


def identify_sky_spectra(spectra, pos, ellipse=None, lmin=650., lmax=700.):

    kt = SedSpec.Spectra(spectra)

    # outer = inner + 3.

    skys = kt.good_positions.tolist()

    a = ellipse[0]*1.25
    b = a * (ellipse[1] / ellipse[0])
    xc = ellipse[2]
    yc = ellipse[3]
    theta = ellipse[4] * (np.pi / 180.)

    all_kix = []
    for the_pos in pos:
        all_kix.append(list(find_positions_ellipse(kt.KT.data, xc, yc, a, b,
                                                   -theta)))
    all_kix = list(itertools.chain(*all_kix))
    kix = list(set(all_kix))
    objs = kt.good_positions[kix]

    for o in objs:
        if o in skys:
            skys.remove(o)

    if len(skys) > 0:
        print "Number of starting pure sky spaxels is %d" % len(skys)
    else:
        print "ERROR: no sky spaxels in this image"

    newspec = [spectra[i] for i in skys]
    kt = SedSpec.Spectra(newspec)

    xs, ys, vs = kt.to_xyv(lmin=lmin, lmax=lmax)
    vmdn = np.median(vs)

    vstd = np.nanstd(vs)

    hi_thresh = vmdn + 1.25 * vstd
    lo_thresh = vmdn - 2.0 * vstd
    print("Median: %6.2f, STD: %6.2f, Hi Thresh: %6.2f, Lo Thresh: %6.2f" %
          (vmdn, vstd, hi_thresh, lo_thresh))

    n_hi_rem = 0
    n_lo_rem = 0
    n_tot = 0

    for s in skys:
        el = spectra[s]
        l, fl = el.get_flambda()

        ok = (l > lmin) & (l <= lmax)

        if np.median(el.spec[ok]) > hi_thresh:
            skys.remove(s)
            n_hi_rem += 1

        if np.median(el.spec[ok]) < lo_thresh:
            skys.remove(s)
            n_lo_rem += 1

        n_tot += 1

    n_tot -= n_hi_rem + n_lo_rem
    print("Removed %d high sky spaxels and %d low sky spaxels leaving %d "
          "remaining spaxels" % (n_hi_rem, n_lo_rem, n_tot))

    return skys


def identify_bgd_spectra(spectra, pos, ellipse=None, expfac=1.):
    kt = SedSpec.Spectra(spectra)

    a = ellipse[0] * expfac
    b = ellipse[1] * expfac
    sky_a = a + 3. * expfac
    sky_b = sky_a * (b/a)
    xc = ellipse[2]
    yc = ellipse[3]
    theta = ellipse[4] * (np.pi / 180.)

    all_kix = []
    for the_pos in pos:
        all_kix.append(list(find_positions_ellipse(kt.KT.data, xc, yc,
                                                   a, b,
                                                   -theta)))
    all_kix = list(itertools.chain(*all_kix))
    kix = list(set(all_kix))
    objs = kt.good_positions[kix]

    all_kix = []
    for the_pos in pos:
        all_kix.append(list(find_positions_ellipse(kt.KT.data, xc, yc,
                                                   sky_a, sky_b,
                                                   -theta)))
    all_kix = list(itertools.chain(*all_kix))
    kix = list(set(all_kix))
    skys = kt.good_positions[kix].tolist()

    for o in objs:
        if o in skys:
            skys.remove(o)

    return skys


def to_image(spectra, meta, outname, posa=None, posb=None, adcpos=None,
             ellipse=None, ellipseb=None, bgd_sub=True,
             lmin=650., lmax=700., cmin=None, cmax=None):
    """ Convert spectra list into image_[outname].pdf

    Args:
        spectra (array of extraction): see Extraction.py
        meta (dict): dictionary of meta-data
        outname (string): output filename
        posa (tuple): x,y position of A aperture (asec)
        posb (tuple): x,y position of B aperture (asec)
        adcpos (tuple): position offsets due to atmospheric dispersion (asec)
        ellipse (tuple): ellipse parameters for A aperture
        ellipseb (tuple): ellipse parameters for B aperture
        bgd_sub (boolean): set to True to subtract median background
        lmin (float): minimum wavelength in nm to sum over
        lmax (float): maximum wavelength in nm to sum over
        cmin (float): cube intensity minimum for scaling
        cmax (float): cube intensity maximum for scaling
    """

    xs = []
    ys = []
    vs = []

    for x in spectra:
        if x.xrange is None:
            continue
        if x.lamcoeff is None:
            continue
        if x.specw is None:
            continue
        ix = np.arange(*x.xrange)
        ll = chebval(ix, x.lamcoeff)
        ok = (ll > lmin) & (ll < lmax)
        if ok.any():
            xs.append(x.X_as)
            ys.append(x.Y_as)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                vs.append(np.median(x.specw[ok]))

    if bgd_sub:
        vs -= np.median(vs)

    # Clean outliers
    vcln = reject_outliers(np.array(vs, dtype=np.float), m=3.)
    vstd = np.nanstd(vcln)
    vmid = np.median(vcln)
    if cmin is None or cmax is None:
        if posb is None:
            vmin = vmid - vstd
            vmax = vmid + 8.*vstd
        else:
            vmin = vmid - 8.*vstd
            vmax = vmid + 8.*vstd
    else:
        vmin = cmin
        vmax = cmax

    pl.clf()
    pl.ylim(-20, 20)
    pl.xlim(-22, 20)
    pl.grid(True)
    if posa is not None:
        pl.axvline(posa[0], color='black', linewidth=.5)
        pl.axhline(posa[1], color='black', linewidth=.5)
    if posb is not None:
        pl.axvline(posb[0], color='black', linewidth=.5)
        pl.axhline(posb[1], color='black', linewidth=.5)
    print "scaling output image between %d and %d" % (vmin, vmax)
    pl.scatter(xs, ys, c=vs, s=50, marker='H', linewidth=0,
               vmin=vmin, vmax=vmax)

    if ellipse is not None:
        xys = get_ellipse_xys(ellipse)
        pl.plot(xys[:, 0], xys[:, 1], 'g.-')

    if ellipseb is not None:
        xys = get_ellipse_xys(ellipseb)
        pl.plot(xys[:, 0], xys[:, 1], 'g.-')

    if adcpos is not None:
        for p in adcpos:
            pl.plot(p[0], p[1], 'rx')

    pl.xlabel("X [as] @ %6.1f nm" % meta['fiducial_wavelength'])
    pl.ylabel("Y [as]")
    tlab = meta['outname']
    if 'airmass' in meta:
        tlab += ", Airmass: %.3f" % meta['airmass']
    pl.title(tlab)
    pl.colorbar()
    pl.savefig("image_%s.pdf" % outname)
    pl.close()
    print "Wrote image_%s.pdf" % outname


def c_to_nm(coefficients, pix, offset=0.):

    t = coefficients[:]
    t[0] += offset
    return chebval(pix, t)


def interp_spectra(all_spectra, six, sign=1., outname=None, plot=False,
                   corrfile=None, dnm=0., onto=None, sky=False, percent=None):
    """Interp spectra onto common grid

    Args:
        all_spectra:
        six:
        sign:
        plot:
        corrfile:
        outname:
        dnm: Offset (usually for flexure) in nm
        onto:
        sky:
        percent:
    """

    l_grid = onto
    s_grid = []
    f_grid = []
    lamcoeff = None
    # for ix,spectrum in enumerate(all_spectra):
    for ix in six:
        spectrum = all_spectra[ix]

        if sky and all_spectra[ix].is_obj:
            continue

        l, s = spectrum.get_counts(the_spec='specw')
        pix = np.arange(*spectrum.xrange)

        # This is wrong: should give preference to lamcoeff according to Nick
        # if spectrum.mdn_coeff is not None: cs = spectrum.mdn_coeff
        # else: cs = spectrum.lamcoeff

        # This is correct: preference to lamcoeff over mdn_coeff
        if spectrum.lamcoeff is not None:
            cs = spectrum.lamcoeff
        else:
            cs = spectrum.mdn_coeff

        # get wavelengths for spectrum
        l = c_to_nm(cs, pix, offset=dnm)

        # skip short spectra (on or near edge of IFU)
        if l.max() - l.min() < 300:
            continue

        # Positive or negative spectra
        pon = sign

        # Check if our wavelength grid is defined,
        if l_grid is None:
            # use the first set of wavelengths and store
            l_grid = l
            fl = s * pon
            lamcoeff = spectrum.lamcoeff
        else:
            # Interpolate onto our wavelength grid and store
            fun = interp1d(l, s*pon, bounds_error=False, fill_value=0)
            fl = fun(l_grid)

        s_grid.append(fl)

        if percent is not None:
            f_grid.append(np.nansum(fl[(l_grid > 500)*(l_grid < 900)]))

    # Are we trimming at a percentile?
    if percent is not None:
        f_lim = np.percentile(f_grid, percent)
        print "Trimming at %.1f %%, flux = %.1f" % (percent, f_lim)
        if f_lim < 0:
            print "WARNING: If A/B pair, sky level has changed between A and B"
            print "         Consider single mode extraction"
        l_grid = onto
        s_grid = []
        newsix = []
        lamcoeff = None
        # for ix,spectrum in enumerate(all_spectra):
        for ix in six:
            spectrum = all_spectra[ix]

            if sky and all_spectra[ix].is_obj:
                continue

            l, s = spectrum.get_counts(the_spec='specw')
            pix = np.arange(*spectrum.xrange)

            # Preference to lamcoeff over mdn_coeff
            if spectrum.lamcoeff is not None:
                cs = spectrum.lamcoeff
            else:
                cs = spectrum.mdn_coeff

            # get wavelengths for spectrum
            l = c_to_nm(cs, pix, offset=dnm)

            # skip short spectra (on or near edge of IFU)
            if l.max() - l.min() < 300:
                continue

            # Positive or negative spectra
            pon = sign

            # Check if our wavelength grid is defined,
            if l_grid is None:
                # use the first set of wavelengths and store
                l_grid = l
                fl = s * pon
                lamcoeff = spectrum.lamcoeff
            else:
                # Interpolate onto our wavelength grid and store
                fun = interp1d(l, s * pon, bounds_error=False, fill_value=0)
                fl = fun(l_grid)

            f_test = np.nansum(fl[(l_grid > 500)*(l_grid < 900)])

            if f_test > f_lim:
                s_grid.append(fl)
                newsix.append(ix)
            else:
                print "rejected - ix: %d, flx: %.1f" % (ix, f_test)
    else:
        newsix = None

    # average of all spectra selected
    medspec = np.nanmean(s_grid, axis=0)

    # Output figures if requested

    # Spectrum
    pl.figure(3)
    pl.clf()
    pl.step(l_grid, medspec)
    yl = pl.ylim()
    pl.xlim(300., 1200.)
    pl.xlabel('Wavelength [nm]')
    pl.ylabel(r'Spectral irradiance[photon/10 m/nm]')
    pl.title("%s Raw Spectrum" % outname.split('.')[0])
    pl.grid(True)
    if outname is not None:
        pl.savefig("spec_%s" % outname)
        print "Wrote spec_%s" % outname
    if plot:
        pl.show()

    # Spaxel stack image
    pl.figure(2)
    pl.clf()
    s_grid = np.array(s_grid)
    pl.imshow(s_grid, vmin=yl[0], vmax=yl[1])
    pl.xlabel('Wavelength bin [pixel]')
    pl.title("%s Spaxels" % outname.split('.')[0])
    pl.colorbar()
    pl.grid(True)
    if outname is not None:
        pl.savefig("allspec_%s" % outname)
        print "Wrote allspec_%s" % outname
    if plot:
        pl.show()

    # Package results
    doc = """Result contains:
        nm [N float]: Wavelength solution
        ph_10m_nm [N float]: Spectral irradiance in units of phot / 10 min / nm
        spectra [? x K float]: List of all the spectra that participated in
            the formation of ph_10m_nm. By interpolating these objects onto
            a ph_10m_nm and taking the mean, you produce ph_10m_nm
        coefficients [3-5 element float]: Chebyshev coefficents that produce
            nm. Can be evaluated with numpy chebval().
        corrected-spec [N float]: ph_10m_nm * Atmospheric correction, if
            available
        doc: This doc string
        """
    result = [{"nm": l_grid, "ph_10m_nm": medspec, "spectra": s_grid,
              "coefficients": lamcoeff, "doc": doc}]

    # Calibrate output if corrfile specified (this is not usually done)
    cc = None

    # Try to load corrfile
    if corrfile is not None:
        try:
            cc = np.load(corrfile)[0]
        except:
            cc = None

    # Apply correction
    if cc is not None:
        corrfun = chebval(l_grid, cc['coeff'])
        corrfun /= np.nanmin(corrfun)
        corrfun = interp1d(cc['nm'], cc['cor'],
                           bounds_error=False, fill_value=np.nan)
        corrfun = corrfun(l_grid)
        result[0]['corrected-spec'] = medspec * corrfun

        # Output corrected spectrum if requested
        pl.figure(4)
        pl.clf()
        pl.step(l_grid, medspec*corrfun)
        pl.ylim(yl[0], yl[1]*20)
        pl.xlabel('Wavelength [nm]')
        pl.ylabel(r'Spectral irradiance[photon/10 m/nm] x Atm correction')
        pl.title("%s Corr Spectrum" % outname.split('.')[0])
        pl.grid(True)
        if outname is not None:
            pl.savefig("corr_spec_%s" % outname)
            print "Wrote corr_spec_%s" % outname
        if plot:
            pl.show()

    pl.figure(2)

    return result, newsix


def imarith(operand1, op, operand2, result, doairmass=False, doexptime=False):
    from pyraf import iraf
    iraf.images()

    pars = iraf.imarith.getParList()
    iraf.imcombine.unlearn()

    try:
        os.remove(result)
    except:
        pass

    print "%s %s %s -> %s" % (operand1, op, operand2, result)
    iraf.imarith(operand1=operand1, op=op, operand2=operand2, result=result)
    iraf.imarith.setParList(pars)
    if doairmass:
        # Adjust FITS header
        with pf.open(operand1) as f:
            am1 = f[0].header['airmass']
        with pf.open(operand2) as f:
            am2 = f[0].header['airmass']

        of = pf.open(result)
        of[0].header['airmass1'] = am1
        of[0].header['airmass2'] = am2
        of.writeto(result, clobber=True)
    if doexptime:
        # Adjust FITS header
        with pf.open(operand1) as f:
            ex1 = f[0].header['exptime']
        with pf.open(operand2) as f:
            ex2 = f[0].header['exptime']

        of = pf.open(result)
        of[0].header['exptime1'] = ex1
        of[0].header['exptime2'] = ex2
        of.writeto(result, clobber=True)


def gunzip(a, b):
    if a.endswith(".gz"):
        os.system("gunzip %s" % a)
        a = os.path.splitext(a)[0]
    if b.endswith(".gz"):
        os.system("gunzip %s" % b)
        b = os.path.splitext(b)[0]

    return a, b


def gzip(a, b):
    if not a.endswith(".gz"):
        os.system("gzip %s" % a)
    if not b.endswith(".gz"):
        os.system("gzip %s" % b)

    return a + ".gz", b + ".gz"


def add(a, b, outname):
    a, b = gunzip(a, b)
    imarith(a, "+", b, outname)
    gzip(a, b)

    return pf.open(outname)


def addcon(a, b, outname):
    a, b = gunzip(a, b)
    imarith(a, "+", b, outname)
    gzip(a, "junk.gz")

    return pf.open(outname)


def subtract(a, b, outname):
    if os.path.exists(outname):
        return pf.open(outname)

    a, b = gunzip(a, b)
    imarith(a, "-", b, outname, doairmass=True)
    a, b = gzip(a, b)

    return pf.open(outname)


def divide(a, b, outname):
    a, b = gunzip(a, b)
    imarith(a, "/", b, outname)
    gzip(a, b)

    return pf.open(outname)


def bgd_level(extractions):
    """Remove background from extractions"""

    levels = []
    for spectrum in extractions:
        if 'spec' in spectrum.__dict__ and spectrum.spec is not None \
          and spectrum.lamcoeff is not None:

            l, fl = spectrum.get_counts(the_spec='specw')

            levels.append(np.median(fl))

    bgd = np.median(levels)
    sd = np.std(levels)
    pl.plot(levels, 'x')
    pl.axhline(bgd)
    pl.axhline(bgd+sd)
    pl.axhline(bgd-sd)
    pl.ylim(-20*sd-bgd, 20*sd-bgd)
    pl.show()


def handle_flat(flfile, fine, outname=None):
    """Loads IFU Flat frame "flfile" and extracts spectra using "fine".

    Args:
        flfile (string): filename of ifu FITS file to extract from.
        fine (string): filename of NumPy file with locations + wavelength soln.
        outname (string): filename to write results to.

    Returns:
        None
    Raises:
        None
    """

    fine, fmeta = np.load(fine)
    if outname is None:
        outname = "%s" % flfile

    if os.path.isfile(outname+".npy"):
        print "Extractions already exist in %s.npy!" % outname
        print "rm %s.npy # if you want to recreate extractions" % outname
    else:
        print "\nCREATING extractions ..."
        spec = pf.open(flfile)

        print "\nExtracting object spectra"
        ex, meta = \
            Wavelength.wavelength_extract(spec, fine,
                                          filename=outname,
                                          flexure_x_corr_nm=0.,
                                          flexure_y_corr_pix=0.)
        meta['airmass'] = spec[0].header['airmass']
        header = {}
        for k, v in spec[0].header.iteritems():
            try:
                header[k] = v
            except:
                pass
        meta['HA'] = header['HA'] if 'HA' in header else header['TEL_HA']
        meta['DEC'] = header['DEC'] if 'DEC' in header else header['TEL_DEC']
        meta['RA'] = header['RA'] if 'RA' in header else header['TEL_RA']
        meta['PRLLTC'] = \
            header['PRLLTC'] if 'PRLLTC' in header else header['TEL_PA']
        meta['EQUINOX'] = \
            header['EQUINOX'] if 'EQUINOX' in header else header['TELEQX']
        if 'UTC' in header:
            meta['utc'] = header['UTC']
        else:
            tjd = Time(header['JD'], format='jd', scale='utc')
            meta['utc'] = tjd.yday
        meta['fiducial_wavelength'] = fmeta['fiducial_wavelength']

        meta['header'] = header

        meta['exptime'] = spec[0].header['exptime']
        np.save(outname, [ex, meta])
        print "Wrote %s.npy" % outname


def handle_std(stdfile, fine, outname=None, standard=None, offset=None,
               flat_corrections=None, lmin=650., lmax=700.,
               refl=0.9, area=18000.):
    """Loads IFU frame "stdfile" and extracts standard star spectra using "fine".

    Args:
        stdfile (string): filename of ifu FITS file to extract from.
        fine (string): filename of NumPy file with locations + wavelength soln
        outname (string): filename to write results to
        standard (string): name of standard star
        offset (2tuple): X (nm)/Y (pix) shift to apply for flexure correction
        flat_corrections (list): A list of FlatCorrection objects for
            correcting the extraction
        lmin (float): lower wavelength limit for image generation
        lmax (float): upper wavelength limit for image generation
        refl (float): Telescope reflectance factor (assuming .90 for P60)
        area (float): Telescope area (assuming 18000. cm^2 for P60)

    Returns:
        The extracted spectrum, a dictionary:
        {'ph_10m_nm': Observed flux in photon / 10 m / nm integrated
        'nm': Wavelength solution in nm
        'N_spax': Total number of spaxels that created ph_10m_nm
        'skyph': Sky flux in photon / 10 m / nanometer / spaxel
        'std-correction': ratio of observed to reference flux
        'std-maxnm': maximum reference wavelength in nm
        'radius_as': Extraction radius in arcsec
        'pos': X/Y extraction location of spectrum in arcsec}

    Raises:
        None
    """

    # Load wavelength/spatial solution
    fine, fmeta = np.load(fine)
    # Set a default outname if needed
    if outname is None:
        outname = "%s" % stdfile
    # Check for flexure offsets
    if offset is not None:
        ff = np.load(offset)
        flexure_x_corr_nm = ff[0]['dXnm']
        flexure_y_corr_pix = ff[0]['dYpix']
        print "Dx %2.1f nm | Dy %2.1f px" % (ff[0]['dXnm'], ff[0]['dYpix'])
        if 'xfwhm' in ff[0]:
            xfwhm = ff[0]['xfwhm']
            yfwhm = ff[0]['yfwhm']
            print "FWHMx %2.1f nm | FWHMy %2.1f px" % (ff[0]['xfwhm'],
                                                       ff[0]['yfwhm'])
        else:
            xfwhm = 0.
            yfwhm = 0.
    else:
        flexure_x_corr_nm = 0.
        flexure_y_corr_pix = 0.
        xfwhm = 0.
        yfwhm = 0.

    # The spaxel extraction already exist, so load them in
    if os.path.isfile(outname+".npy"):
        print "USING extractions in %s.npy!" % outname
        print "rm %s.npy # if you want to recreate extractions" % outname
        ex, meta = np.load(outname+".npy")
        e_var, meta_var = np.load("var_" + outname + ".npy")
    # No extractions yet, so generate them
    else:
        print "\nCREATING extractions ..."
        # Load spectrum image
        spec = pf.open(stdfile)
        # Set up variance image by adding read noise to Poisson noise
        adcspeed = spec[0].header["ADCSPEED"]
        if adcspeed == 2:
            read_var = 22*22
        else:
            read_var = 5*5
        # Variance image is input image plus read variance
        var = addcon(stdfile, str(read_var), "var_" + outname + ".fits")
        # Extract each object spaxel
        print "\nExtracting object spectra"
        ex, meta = \
            Wavelength.wavelength_extract(spec, fine,
                                          filename=outname,
                                          flexure_x_corr_nm=flexure_x_corr_nm,
                                          flexure_y_corr_pix=flexure_y_corr_pix,
                                          flat_corrections=flat_corrections)
        # Extract metadata
        meta['airmass'] = spec[0].header['airmass']
        header = {}
        for k, v in spec[0].header.iteritems():
            try:
                header[k] = v
            except:
                pass
        meta['HA'] = header['HA'] if 'HA' in header else header['TEL_HA']
        meta['DEC'] = header['DEC'] if 'DEC' in header else header['TEL_DEC']
        meta['RA'] = header['RA'] if 'RA' in header else header['TEL_RA']
        meta['PRLLTC'] = \
            header['PRLLTC'] if 'PRLLTC' in header else header['TEL_PA']
        meta['EQUINOX'] = \
            header['EQUINOX'] if 'EQUINOX' in header else header['TELEQX']
        if 'UTC' in header:
            meta['utc'] = header['UTC']
        else:
            tjd = Time(header['JD'], format='jd', scale='utc')
            meta['utc'] = tjd.yday
        meta['fiducial_wavelength'] = fmeta['fiducial_wavelength']
        meta['header'] = header
        meta['exptime'] = spec[0].header['exptime']
        # Save object extraction
        np.save(outname, [ex, meta])
        print "Wrote %s.npy" % outname
        # Extract each variance spaxel
        print "\nExtracting variance spectra"
        e_var, meta_var = \
            Wavelength.wavelength_extract(var, fine,
                                          filename=outname,
                                          flexure_x_corr_nm=flexure_x_corr_nm,
                                          flexure_y_corr_pix=flexure_y_corr_pix,
                                          flat_corrections=flat_corrections)
        # Save variance extraction
        np.save("var_" + outname, [e_var, meta_var])
        print "Wrote var_%s.npy" % outname

    # Get the object name of record
    objname = meta['header']['OBJECT'].split()[0]

    # Automatic extraction using Gaussian fit for Standard Stars
    sixa, posa, adcpos, ellipse, status = \
        identify_spectra_gauss_fit(ex,
                                   prlltc=Angle(meta['PRLLTC'], unit='deg'),
                                   lmin=lmin, lmax=lmax,
                                   airmass=meta['airmass'])

    if status > 0:
        quality = 4  # something went wrong
    else:
        quality = 0  # good fit

    radius_used = ellipse[0] * 0.5

    # Mark object spaxels
    for ix in sixa:
        ex[ix].is_obj = True
    # Use all sky spaxels in image for Standard Stars
    kixa = identify_sky_spectra(ex, posa, ellipse=ellipse)
    # Mark sky spaxels
    for ix in kixa:
        ex[ix].is_sky = True

    # Make an image of the spaxels for the record
    to_image(ex, meta, outname, posa=posa, adcpos=adcpos, ellipse=ellipse,
             bgd_sub=False, lmin=lmin, lmax=lmax)
    # get the mean spectrum over the selected spaxels
    resa, nspxa = interp_spectra(ex, sixa, outname=outname+".pdf")
    skya, nspxak = interp_spectra(ex, kixa, outname=outname+"_sky.pdf",
                                  sky=True)
    vara, nspxav = interp_spectra(e_var, sixa, outname=outname+"_var.pdf")

    # Plot out the X/Y positions of the selected spaxels
    xsa = []
    ysa = []
    xsk = []
    ysk = []
    for ix in sixa:
        xsa.append(ex[ix].X_as)
        ysa.append(ex[ix].Y_as)
    for ix in kixa:
        if not ex[ix].is_obj:
            xsk.append(ex[ix].X_as)
            ysk.append(ex[ix].Y_as)

    pl.figure()
    pl.clf()
    pl.ylim(-20, 20)
    pl.xlim(-22, 20)
    pl.grid(True)

    if posa is not None:
        pl.axvline(posa[0], color='black', linewidth=.5)
        pl.axhline(posa[1], color='black', linewidth=.5)
    if ellipse is not None:
        xys = get_ellipse_xys(ellipse)
        pl.plot(xys[:, 0], xys[:, 1], 'g.-')

    pl.xlabel("X [as] @ %6.1f nm" % meta['fiducial_wavelength'])
    pl.ylabel("Y [as]")
    pl.scatter(xsa, ysa, color='red', marker='H', s=50, linewidth=0)
    pl.scatter(xsk, ysk, color='green', marker='H', s=50, linewidth=0)
    tlab = "%d selected spaxels for %s" % (len(xsa), objname)
    if 'airmass' in meta:
        tlab += "\nAirmass: %.3f" % meta['airmass']
    pl.title(tlab)
    pl.savefig("XYs_%s.pdf" % outname)
    pl.close()
    print "Wrote XYs_%s.pdf" % outname
    # / End Plot

    # Define our standard wavelength grid
    ll = Wavelength.fiducial_spectrum()
    # Resample sky onto standard wavelength grid
    sky_a = interp1d(skya[0]['nm'], skya[0]['ph_10m_nm'], bounds_error=False)
    sky = sky_a(ll)
    # Resample variance onto standard wavelength grid
    var_a = interp1d(vara[0]['nm'], vara[0]['ph_10m_nm'], bounds_error=False)
    varspec = var_a(ll)
    # Copy and resample object spectrum onto standard wavelength grid
    res = [{"doc": resa[0]["doc"], "ph_10m_nm": np.copy(resa[0]["ph_10m_nm"]),
            "spectra": np.copy(resa[0]["spectra"]),
            "coefficients": np.copy(resa[0]["coefficients"]),
            "nm": np.copy(resa[0]["ph_10m_nm"])}]
    res[0]['nm'] = np.copy(ll)
    f1 = interp1d(resa[0]['nm'], resa[0]['ph_10m_nm'], bounds_error=False)
    # Calculate airmass correction
    airmass = meta['airmass']
    extcorr = 10**(Atm.ext(ll*10) * airmass/2.5)
    print "Median airmass corr: %.4f" % np.median(extcorr)
    # Calculate output corrected spectrum
    # Account for sky, airmass and aperture
    res[0]['ph_10m_nm'] = (f1(ll)-sky_a(ll)) * extcorr * len(sixa)

    # Process standard star objects
    print "STANDARD"
    # Extract reference data
    wav = standard[:, 0]
    flux = standard[:, 1]
    # Flux in photons/s/cm^2/A
    flpho = 5.0341125e-9 * flux * wav
    # convert wav to nm
    wav /= 10.
    fun = interp1d(wav, flpho, bounds_error=False, fill_value=np.nan)
    # Effective area
    earea = (res[0]['ph_10m_nm'] / 600.) / (fun(res[0]['nm']) * 10.)
    # Efficiency assuming given reflectance (refl) and area in cm^2
    eff = earea * refl / area
    # Calculate/Interpolate correction onto object wavelengths
    fun = interp1d(wav, flux, bounds_error=False, fill_value=np.nan)
    # Divide reference spectrum by observed to get correction
    correction0 = fun(res[0]['nm']) / res[0]['ph_10m_nm']
    # Filter for resolution
    flxf = filters.gaussian_filter(flux, 19.)
    # Calculate/Interpolate filtered correction
    fun = interp1d(wav, flxf, bounds_error=False, fill_value=np.nan)
    # Divide reference spectrum by observed to get correction
    correction = fun(res[0]['nm']) / res[0]['ph_10m_nm']
    # Use unfiltered for H-beta region
    roi = (res[0]['nm'] > 470.) & (res[0]['nm'] < 600.)
    correction[roi] = correction0[roi]
    # Store correction and max calibrated wavelength
    res[0]['std-correction'] = correction
    res[0]['std-maxnm'] = np.max(wav)
    res[0]['reflectance'] = refl
    res[0]['area'] = area
    res[0]['ea'] = earea
    res[0]['efficiency'] = eff

    # Store flexure data
    res[0]['dXnm'] = flexure_x_corr_nm
    res[0]['dYpix'] = flexure_y_corr_pix
    res[0]['xfwhm'] = xfwhm
    res[0]['yfwhm'] = yfwhm

    # Store new metadata
    res[0]['exptime'] = meta['exptime']
    res[0]['Extinction Correction'] = 'Applied using Hayes & Latham'
    res[0]['extinction_corr'] = extcorr
    res[0]['skyph'] = sky * len(sixa)
    res[0]['skynm'] = ll
    res[0]['var'] = varspec
    res[0]['radius_as'] = radius_used
    res[0]['position'] = posa
    res[0]['N_spax'] = len(sixa)
    res[0]['meta'] = meta
    res[0]['object_spaxel_ids'] = sixa
    res[0]['sky_spaxel_ids'] = kixa
    res[0]['sky_spectra'] = skya[0]['spectra']
    res[0]['sky_subtraction'] = True
    res[0]['quality'] = quality
    # Calculate wavelength offsets
    coef = chebfit(np.arange(len(ll)), ll, 4)
    xs = np.arange(len(ll)+1)
    newll = chebval(xs, coef)
    # Store offsets
    res[0]['dlam'] = np.diff(newll)
    # Save the final spectrum
    np.save("sp_" + outname, res)
    print "Wrote sp_"+outname+".npy"
    # Save std star ellipse if all went well
    if quality == 0:
        np.save("ell_" + outname, ellipse)
        print "Wrote ell_%s.npy" % outname


def handle_single(imfile, fine, outname=None, offset=None,
                  radius=2., flat_corrections=None, nosky=False,
                  lmin=650., lmax=700.,
                  specExtract=False, autoExtract=False, interact=False):
    """Loads IFU frame "imfile" and extracts spectra using "fine".

    Args:
        imfile (string): filename of ifu FITS file to extract from.
        fine (string): filename of NumPy file with locations + wavelength
            soln
        outname (string): filename to write results to
        offset (2tuple): X (nm)/Y (pix) shift to apply for flexure correction
        radius (float): Extraction radius in arcsecond
        flat_corrections (list): A list of FlatCorrection objects for
            correcting the extraction
        nosky (Boolean): if True don't subtract sky, merely sum in aperture
        lmin (float): lower wavelength limit for image generation
        lmax (float): upper wavelength limit for image generation
        specExtract (Boolean): perform extraction to a spectrum?
        autoExtract (Boolean): automatically find source?
        interact (Boolean): Interactively set the quality of the extraction?

    Returns:
        The extracted spectrum, a dictionary:
        {'ph_10m_nm': Flux in photon / 10 m / nanometer integrated
        'nm': Wavelength solution in nm
        'N_spax': Total number of spaxels that created ph_10m_nm
        'skyph': Sky flux in photon / 10 m / nanometer / spaxel
        'radius_as': Extraction radius in arcsec
        'pos': X/Y extraction location of spectrum in arcsec}

    Raises:
        None
    """

    # Load wavelength/spatial solution
    fine, fmeta = np.load(fine)
    # Set a default outname if needed
    if outname is None:
        outname = "%s" % imfile
    # Check for flexure offsets
    if offset is not None:
        ff = np.load(offset)
        flexure_x_corr_nm = ff[0]['dXnm']
        flexure_y_corr_pix = ff[0]['dYpix']
        print "Dx %2.1f nm | Dy %2.1f px" % (ff[0]['dXnm'], ff[0]['dYpix'])
        if 'xfwhm' in ff[0]:
            xfwhm = ff[0]['xfwhm']
            yfwhm = ff[0]['yfwhm']
            print "FWHMx %2.1f nm | FWHMy %2.1f px" % (ff[0]['xfwhm'],
                                                       ff[0]['yfwhm'])
        else:
            xfwhm = 0.
            yfwhm = 0.
    else:
        flexure_x_corr_nm = 0
        flexure_y_corr_pix = 0
        xfwhm = 0.
        yfwhm = 0.

    # The spaxel extraction already exist, so load them in
    if os.path.isfile(outname+".npy"):
        print "USING extractions in %s.npy!" % outname
        print "rm %s.npy # if you want to recreate extractions" % outname
        ex, meta = np.load(outname+".npy")
        e_var, meta_var = np.load("var_" + outname + ".npy")
    # No extractions yet, so generate them
    else:
        print "\nCREATING extractions ..."
        # Load spectrum image
        spec = pf.open(imfile)
        # Set up variance image by adding read noise to Poisson noise
        adcspeed = spec[0].header["ADCSPEED"]
        if adcspeed == 2:
            read_var = 22*22
        else:
            read_var = 5*5
        # Variance image is input image plus read variance
        var = addcon(imfile, str(read_var), "var_" + outname + ".fits")
        # Extract each object spaxel
        print "\nExtracting object spectra"
        ex, meta = \
            Wavelength.wavelength_extract(spec, fine,
                                          filename=outname,
                                          flexure_x_corr_nm=flexure_x_corr_nm,
                                          flexure_y_corr_pix=flexure_y_corr_pix,
                                          flat_corrections=flat_corrections)
        # Extract metadata
        meta['airmass'] = spec[0].header['airmass']
        header = {}
        for k, v in spec[0].header.iteritems():
            try:
                header[k] = v
            except:
                pass
        meta['HA'] = header['HA'] if 'HA' in header else header['TEL_HA']
        meta['DEC'] = header['DEC'] if 'DEC' in header else header['TEL_DEC']
        meta['RA'] = header['RA'] if 'RA' in header else header['TEL_RA']
        meta['PRLLTC'] = \
            header['PRLLTC'] if 'PRLLTC' in header else header['TEL_PA']
        meta['EQUINOX'] = \
            header['EQUINOX'] if 'EQUINOX' in header else header['TELEQX']
        if 'UTC' in header:
            meta['utc'] = header['UTC']
        else:
            tjd = Time(header['JD'], format='jd', scale='utc')
            meta['utc'] = tjd.yday
        meta['fiducial_wavelength'] = fmeta['fiducial_wavelength']
        meta['header'] = header
        meta['exptime'] = spec[0].header['exptime']
        # Save object extraction
        np.save(outname, [ex, meta])
        print "Wrote %s.npy" % outname
        # Extract each variance spaxel
        print "\nExtracting variance spectra"
        e_var, meta_var = \
            Wavelength.wavelength_extract(var, fine,
                                          filename=outname,
                                          flexure_x_corr_nm=flexure_x_corr_nm,
                                          flexure_y_corr_pix=flexure_y_corr_pix,
                                          flat_corrections=flat_corrections)
        # Save variance extraction
        np.save("var_" + outname, [e_var, meta_var])
        print "Wrote var_%s.npy" % outname

    if specExtract:
        # Get the object name of record
        objname = meta['header']['OBJECT'].split()[0]

        if autoExtract:
            # Automatic extraction using Gaussian fit for Standard Stars
            sixa, posa, adcpos, ellipse, status = \
                identify_spectra_gauss_fit(ex,
                                           prlltc=Angle(meta['PRLLTC'],
                                                        unit='deg'),
                                           lmin=lmin, lmax=lmax,
                                           airmass=meta['airmass'],
                                           sigfac=2.0)
            radius_used = ellipse[0] * 0.5
            if status > 0:
                quality = 4  # something went wrong
            else:
                quality = 1  # by definition

            # Use all sky spaxels in image
            kixa = identify_sky_spectra(ex, posa, ellipse=ellipse)

        else:
            message = "\nMark positive (red) target"

            # A single-frame Science Object
            sixa, posa, adcpos, ellipse, stats = \
                identify_spectra_gui(ex, radius=radius,
                                     prlltc=Angle(meta['PRLLTC'], unit='deg'),
                                     scaled=False, bgd_sub=False,
                                     lmin=lmin, lmax=lmax,
                                     objname=objname, airmass=meta['airmass'],
                                     nosky=nosky,
                                     message=message)
            radius_used = ellipse[0]

            if interact:

                # Get quality of observation
                print "Enter quality of observation:"
                print "1 - good       (no problems)"
                print "2 - acceptable (minor problem)"
                print "3 - poor       (major problem)"
                print "4 - no object visible"
                q = 'x'
                quality = -1
                prom = ": "
                while not q.isdigit() or quality < 1 or quality > 4:
                    q = raw_input(prom)
                    if q.isdigit():
                        quality = int(q)
                        if quality < 1 or quality > 4:
                            prom = "Try again: "
                    else:
                        prom = "Try again: "
                print "Quality = %d, now making outputs..." % quality
            else:
                quality = 0
                print "Now making outputs..."

            # Use an annulus for sky spaxels for Science Objects
            kixa = identify_bgd_spectra(ex, posa, ellipse=ellipse, expfac=1.5)

        for ix in sixa:
            ex[ix].is_obj = True
        for ix in kixa:
            ex[ix].is_sky = True

        # Make an image of the spaxels for the record
        to_image(ex, meta, outname, posa=posa, adcpos=adcpos, ellipse=ellipse,
                 bgd_sub=False, lmin=lmin, lmax=lmax)
        # get the mean spectrum over the selected spaxels
        resa, nsxa = interp_spectra(ex, sixa, outname=outname+".pdf",
                                    percent=30.)
        skya, nsxak = interp_spectra(ex, kixa, outname=outname+"_sky.pdf",
                                     sky=True)
        vara, nsxav = interp_spectra(e_var, nsxa, outname=outname+"_var.pdf")
        # Plot out the X/Y positions of the selected spaxels
        xsa = []
        ysa = []
        xsk = []
        ysk = []
        for ix in nsxa:
            xsa.append(ex[ix].X_as)
            ysa.append(ex[ix].Y_as)
        for ix in kixa:
            if not ex[ix].is_obj:
                xsk.append(ex[ix].X_as)
                ysk.append(ex[ix].Y_as)

        pl.figure()
        pl.clf()
        pl.ylim(-20, 20)
        pl.xlim(-22, 20)
        pl.grid(True)

        if posa is not None:
            pl.axvline(posa[0], color='black', linewidth=.5)
            pl.axhline(posa[1], color='black', linewidth=.5)
        if ellipse is not None:
            xys = get_ellipse_xys(ellipse)
            pl.plot(xys[:, 0], xys[:, 1], 'g.-')

        pl.xlabel("X [as] @ %6.1f nm" % meta['fiducial_wavelength'])
        pl.ylabel("Y [as]")
        pl.scatter(xsa, ysa, color='red', marker='H', s=50, linewidth=0)
        pl.scatter(xsk, ysk, color='green', marker='H', s=50, linewidth=0)
        tlab = "%d selected spaxels for %s" % (len(xsa), objname)
        if 'airmass' in meta:
            tlab += "\nAirmass: %.3f" % meta['airmass']
        if 1 <= quality <= 4:
            tlab += ", Qual: %d" % quality
        pl.title(tlab)
        pl.savefig("XYs_%s.pdf" % outname)
        print "Wrote XYs_%s.pdf" % outname
        pl.close()
        # / End Plot

        # Define our standard wavelength grid
        ll = Wavelength.fiducial_spectrum()
        # Resample sky onto standard wavelength grid
        sky_a = interp1d(skya[0]['nm'], skya[0]['ph_10m_nm'],
                         bounds_error=False)
        sky = sky_a(ll)
        # Resample variance onto standard wavelength grid
        var_a = interp1d(vara[0]['nm'], vara[0]['ph_10m_nm'],
                         bounds_error=False)
        varspec = var_a(ll)
        # Copy and resample object spectrum onto standard wavelength grid
        res = [{"doc": resa[0]["doc"],
                "ph_10m_nm": np.copy(resa[0]["ph_10m_nm"]),
                "spectra": np.copy(resa[0]["spectra"]),
                "coefficients": np.copy(resa[0]["coefficients"]),
                "nm": np.copy(resa[0]["ph_10m_nm"])}]
        res[0]['nm'] = np.copy(ll)
        f1 = interp1d(resa[0]['nm'], resa[0]['ph_10m_nm'], bounds_error=False)
        # Calculate airmass correction
        airmass = meta['airmass']
        extcorr = 10**(Atm.ext(ll*10) * airmass/2.5)
        print "Median airmass corr: %.4f" % np.median(extcorr)
        # Calculate output corrected spectrum
        if stats['nosky']:
            print "Sky subtraction off"
            # Account for airmass and aperture
            res[0]['ph_10m_nm'] = f1(ll) * extcorr * len(nsxa)
        else:
            print "Sky subtraction on"
            # Account for sky, airmass and aperture
            res[0]['ph_10m_nm'] = (f1(ll)-sky_a(ll)) * extcorr * len(nsxa)

        # Store flexure data
        res[0]['dXnm'] = flexure_x_corr_nm
        res[0]['dYpix'] = flexure_y_corr_pix
        res[0]['xfwhm'] = xfwhm
        res[0]['yfwhm'] = yfwhm

        # Store new metadata
        res[0]['exptime'] = meta['exptime']
        res[0]['Extinction Correction'] = 'Applied using Hayes & Latham'
        res[0]['extinction_corr'] = extcorr
        res[0]['skyph'] = sky * len(nsxa)
        res[0]['skynm'] = ll
        res[0]['var'] = varspec
        res[0]['radius_as'] = radius_used
        res[0]['position'] = posa
        res[0]['N_spax'] = len(nsxa)
        res[0]['meta'] = meta
        res[0]['object_spaxel_ids'] = nsxa
        res[0]['sky_spaxel_ids'] = kixa
        res[0]['sky_spectra'] = skya[0]['spectra']
        res[0]['sky_subtraction'] = False if stats['nosky'] else True
        res[0]['quality'] = quality
        # Calculate wavelength offsets
        coef = chebfit(np.arange(len(ll)), ll, 4)
        xs = np.arange(len(ll)+1)
        newll = chebval(xs, coef)
        # Store offsets
        res[0]['dlam'] = np.diff(newll)
        # Save the final spectrum
        np.save("sp_" + outname, res)
        print "Wrote sp_"+outname+".npy"


def handle_dual(afile, bfile, fine, outname=None, offset=None, radius=2.,
                flat_corrections=None, nosky=False, lmin=650., lmax=700.,
                specExtract=False, interact=False):
    """Loads IFU frame "afile" and "bfile" and extract A-B spectra using "fine".

    Args:
        afile (string): filename of ifu FITS file to extract from.
        bfile (string): filename of ifu FITS file to extract from.
        fine (string): filename of NumPy file with locations + wavelength soln.
        outname (string): filename to write results to
        offset (2tuple): X (nm)/Y (pix) shift to apply for flexure correction
        radius (float): Extraction radius in arcsecond
        flat_corrections (list): A list of FlatCorrection objects for
            correcting the extraction
        nosky (Boolean): if True don't subtract sky, merely sum in aperture
        lmin (float): lower wavelength limit for image generation
        lmax (float): upper wavelength limit for image generation
        specExtract (Boolean): perform extraction to a spectrum?
        interact (Boolean): Interactively set the quality of the extraction?

    Returns:
        The extracted spectrum, a dictionary:
        {'ph_10m_nm': Flux in photon / 10 m / nanometer integrated
        'var'
        'nm': Wavelength solution in nm
        'N_spaxA': Total number of "A" spaxels
        'N_spaxB': Total number of "B" spaxels
        'skyph': Sky flux in photon / 10 m / nanometer / spaxel
        'radius_as': Extraction radius in arcsec
        'pos': X/Y extraction location of spectrum in arcsec}

    Raises:
        None
    """

    fine, fmeta = np.load(fine)
    if outname is None:
        outname = "%sm%s" % (afile, bfile)

    if offset is not None:
        ff = np.load(offset)
        flexure_x_corr_nm = ff[0]['dXnm']
        flexure_y_corr_pix = -ff[0]['dYpix']
        print "Dx %2.1f nm | Dy %2.1f px" % (ff[0]['dXnm'], ff[0]['dYpix'])
        if 'xfwhm' in ff[0]:
            xfwhm = ff[0]['xfwhm']
            yfwhm = ff[0]['yfwhm']
            print "FWHMx %2.1f nm | FWHMy %2.1f px" % (ff[0]['xfwhm'],
                                                       ff[0]['yfwhm'])
        else:
            xfwhm = 0.
            yfwhm = 0.
    else:
        flexure_x_corr_nm = 0
        flexure_y_corr_pix = 0
        xfwhm = 0.
        yfwhm = 0.

    if os.path.isfile(outname + ".npy"):
        print "USING extractions in %s!" % outname
        print "rm %s.npy # if you want to recreate extractions" % outname
        ex, meta = np.load(outname + ".npy")
        ex_var, meta_var = np.load("var_" + outname + ".npy")
        header = meta['header']
    else:
        print "\nCREATING extractions ..."
        diff = subtract(afile, bfile, outname + ".fits")
        add(afile, bfile, "tmpvar_" + outname + ".fits")

        adcspeed = diff[0].header["ADCSPEED"]
        if adcspeed == 2:
            read_var = 22*22
        else:
            read_var = 5*5

        var = addcon("tmpvar_" + outname + ".fits", str(read_var),
                     "var_" + outname + ".fits")
        os.remove("tmpvar_" + outname + ".fits.gz")

        print "\nExtracting object spectra"
        ex, meta = \
            Wavelength.wavelength_extract(diff, fine,
                                          filename=outname,
                                          flexure_x_corr_nm=flexure_x_corr_nm,
                                          flexure_y_corr_pix=flexure_y_corr_pix,
                                          flat_corrections=flat_corrections)
        meta['airmass1'] = diff[0].header['airmass1']
        meta['airmass2'] = diff[0].header['airmass2']
        meta['airmass'] = diff[0].header['airmass']
        header = {}
        for k, v in diff[0].header.iteritems():
            try:
                header[k] = v
            except:
                pass
        meta['HA'] = header['HA'] if 'HA' in header else header['TEL_HA']
        meta['DEC'] = header['DEC'] if 'DEC' in header else header['TEL_DEC']
        meta['RA'] = header['RA'] if 'RA' in header else header['TEL_RA']
        meta['PRLLTC'] = \
            header['PRLLTC'] if 'PRLLTC' in header else header['TEL_PA']
        meta['EQUINOX'] = \
            header['EQUINOX'] if 'EQUINOX' in header else header['TELEQX']
        if 'UTC' in header:
            meta['utc'] = header['UTC']
        else:
            tjd = Time(header['JD'], format='jd', scale='utc')
            meta['utc'] = tjd.yday
        meta['fiducial_wavelength'] = fmeta['fiducial_wavelength']

        meta['header'] = header

        meta['exptime'] = diff[0].header['exptime']
        np.save(outname, [ex, meta])
        print "Wrote %s.npy" % outname

        print "\nExtracting variance spectra"
        ex_var, meta_var = \
            Wavelength.wavelength_extract(var, fine,
                                          filename=outname,
                                          flexure_x_corr_nm=flexure_x_corr_nm,
                                          flexure_y_corr_pix=flexure_y_corr_pix,
                                          flat_corrections=flat_corrections)

        np.save("var_" + outname, [ex_var, meta_var])
        print "Wrote var_%s.npy" % outname

    if specExtract:
        objname = header['OBJECT'].split()[0]

        message = "\nMark positive (red) target first"

        sixa, posa, adc_a, ellipse, stats = \
            identify_spectra_gui(ex, radius=radius,
                                 prlltc=Angle(meta['PRLLTC'], unit='deg'),
                                 scaled=False,
                                 lmin=lmin, lmax=lmax,
                                 objname=objname, airmass=meta['airmass'],
                                 nosky=nosky,
                                 message=message)
        radius_used_a = ellipse[0]
        for ix in sixa:
            ex[ix].is_obj = True

        message = "\nMark negative (blue) target next"

        sixb, posb, adc_b, ellipseb, stats = \
            identify_spectra_gui(ex, radius=radius_used_a,
                                 prlltc=Angle(meta['PRLLTC'], unit='deg'),
                                 scaled=stats["scaled"],
                                 lmin=stats["lmin"], lmax=stats["lmax"],
                                 cmin=stats["cmin"], cmax=stats["cmax"],
                                 objname=objname, airmass=meta['airmass'],
                                 nosky=stats["nosky"],
                                 message=message)
        for ix in sixb:
            ex[ix].is_obj = True

        if interact:

            # Get quality of observation
            print "Enter quality of observation:"
            print "1 - good       (no problems)"
            print "2 - acceptable (minor problem)"
            print "3 - poor       (major problem)"
            print "4 - no object visible"
            q = 'x'
            quality = -1
            prom = ": "
            while not q.isdigit() or quality < 1 or quality > 4:
                q = raw_input(prom)
                if q.isdigit():
                    quality = int(q)
                    if quality < 1 or quality > 4:
                        prom = "Try again: "
                else:
                    prom = "Try again: "
            print "Quality = %d, now making outputs..." % quality
        else:
            quality = 0
            print "Now making outputs..."

        to_image(ex, meta, outname, posa=posa, posb=posb, adcpos=adc_a,
                 ellipse=ellipse, ellipseb=ellipseb,
                 lmin=lmin, lmax=lmax,
                 cmin=stats["cmin"],
                 cmax=stats["cmax"])

        kixa = identify_bgd_spectra(ex, posa, ellipse=ellipse)
        for ix in kixa:
            ex[ix].is_sky = True
            if ex[ix].is_obj:
                kixa.remove(ix)
        kixb = identify_bgd_spectra(ex, posb, ellipse=ellipseb)
        for ix in kixb:
            ex[ix].is_sky = True
            if ex[ix].is_obj:
                kixb.remove(ix)

        resa, nsxA = interp_spectra(ex, sixa, sign=1,
                                    outname=outname+"_A.pdf", percent=30.)
        resb, nsxB = interp_spectra(ex, sixb, sign=-1,
                                    outname=outname+"_B.pdf", percent=30.)
        skya, nsxAk = interp_spectra(ex, kixa, sign=1,
                                     outname=outname+"_skyA.pdf", sky=True)
        skyb, nsxBk = interp_spectra(ex, kixb, sign=-1,
                                     outname=outname+"_skyB.pdf", sky=True)
        vara, nsxAv = interp_spectra(ex_var, nsxA, sign=1,
                              outname=outname+"_A_var.pdf")
        varb, nsxBv = interp_spectra(ex_var, nsxB, sign=1,
                              outname=outname+"_B_var.pdf")

        # Plot out the X/Y selected spectra
        xsa = []
        ysa = []
        xsb = []
        ysb = []
        xka = []
        yka = []
        xkb = []
        ykb = []
        for ix in nsxA:
            xsa.append(ex[ix].X_as)
            ysa.append(ex[ix].Y_as)
        for ix in nsxB:
            xsb.append(ex[ix].X_as)
            ysb.append(ex[ix].Y_as)
        for ix in kixa:
            if not ex[ix].is_obj:
                xka.append(ex[ix].X_as)
                yka.append(ex[ix].Y_as)
        for ix in kixb:
            if not ex[ix].is_obj:
                xkb.append(ex[ix].X_as)
                ykb.append(ex[ix].Y_as)

        pl.figure()
        pl.clf()
        pl.ylim(-20, 20)
        pl.xlim(-22, 20)
        pl.grid(True)

        if posa is not None:
            pl.axvline(posa[0], color='black', linewidth=.5)
            pl.axhline(posa[1], color='black', linewidth=.5)
        if posb is not None:
            pl.axvline(posb[0], color='black', linewidth=.5)
            pl.axhline(posb[1], color='black', linewidth=.5)

        if ellipse is not None:
            xys = get_ellipse_xys(ellipse)
            pl.plot(xys[:, 0], xys[:, 1], 'g.-')

        if ellipseb is not None:
            xys = get_ellipse_xys(ellipseb)
            pl.plot(xys[:, 0], xys[:, 1], 'g.-')

        pl.xlabel("X [as] @ %6.1f nm" % meta['fiducial_wavelength'])
        pl.ylabel("Y [as]")
        tlab = "%d selected spaxels for %s" % ((len(nsxA) + len(nsxB)),
                                               meta['outname'])
        air = (meta['airmass1'] + meta['airmass2']) / 2.
        tlab += ", Airmass: %.3f" % air
        if 1 <= quality <= 4:
            tlab += ", Qual: %d" % quality
        pl.title(tlab)
        pl.scatter(xsa, ysa, color='red', marker='H', s=50, linewidth=0)
        pl.scatter(xsb, ysb, color='blue', marker='H', s=50, linewidth=0)
        pl.scatter(xka, yka, color='green', marker='H', s=50, linewidth=0)
        pl.scatter(xkb, ykb, color='green', marker='H', s=50, linewidth=0)
        pl.savefig("XYs_%s.pdf" % outname)
        print "Wrote XYs_%s.pdf" % outname
        pl.close()
        # / End Plot

        np.save("sp_A_" + outname, resa)
        np.save("sp_B_" + outname, resb)
        np.save("var_A_" + outname, vara)
        np.save("var_B_" + outname, varb)
        print("Wrote sp_A_%s.npy, sp_B_%s.npy, var_A_%s.npy, var_B_%s.npy" %
              (outname, outname, outname, outname))

        ll = Wavelength.fiducial_spectrum()
        sky_a = interp1d(skya[0]['nm'], skya[0]['ph_10m_nm'],
                         bounds_error=False)
        sky_b = interp1d(skyb[0]['nm'], skyb[0]['ph_10m_nm'],
                         bounds_error=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            sky = np.nanmean([sky_a(ll), sky_b(ll)], axis=0)

        var_a = interp1d(vara[0]['nm'], vara[0]['ph_10m_nm'],
                         bounds_error=False)
        var_b = interp1d(varb[0]['nm'], varb[0]['ph_10m_nm'],
                         bounds_error=False)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            varspec = np.nanmean([var_a(ll), var_b(ll)], axis=0) * \
                      (len(nsxA) + len(nsxB))

        res = [{"doc": resa[0]["doc"],
                "ph_10m_nm": np.copy(resa[0]["ph_10m_nm"]),
                "nm": np.copy(resa[0]["ph_10m_nm"])}]
        res[0]['nm'] = np.copy(ll)
        f1 = interp1d(resa[0]['nm'], resa[0]['ph_10m_nm'], bounds_error=False)
        f2 = interp1d(resb[0]['nm'], resb[0]['ph_10m_nm'], bounds_error=False)

        airmassa = meta['airmass1']
        airmassb = meta['airmass2']

        extcorra = 10**(Atm.ext(ll*10)*airmassa/2.5)
        extcorrb = 10**(Atm.ext(ll*10)*airmassb/2.5)
        print("Median airmass corrs A: %.4f, B: %.4f" %
              (np.median(extcorra), np.median(extcorrb)))
        # If requested merely sum in aperture, otherwise subtract sky
        if stats['nosky']:
            print "Sky subtraction off"
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=FutureWarning)
                res[0]['ph_10m_nm'] = \
                    np.nansum([f1(ll) * extcorra, f2(ll) * extcorrb],
                              axis=0) * (len(nsxA) + len(nsxB))
        else:
            print "Sky subtraction on"
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=FutureWarning)
                res[0]['ph_10m_nm'] = \
                    np.nansum([(f1(ll)-sky_a(ll)) * extcorra,
                               (f2(ll)-sky_b(ll)) * extcorrb], axis=0) * \
                             (len(nsxA) + len(nsxB))

        # Store flexure data
        res[0]['dXnm'] = flexure_x_corr_nm
        res[0]['dYpix'] = flexure_y_corr_pix
        res[0]['xfwhm'] = xfwhm
        res[0]['yfwhm'] = yfwhm

        # Store new metatdata
        res[0]['exptime'] = meta['exptime']
        res[0]['Extinction Correction'] = 'Applied using Hayes & Latham'
        res[0]['extinction_corr_A'] = extcorra
        res[0]['extinction_corr_B'] = extcorrb
        res[0]['skyph'] = sky * (len(nsxA) + len(nsxB))
        res[0]['var'] = varspec
        res[0]['radius_as'] = radius_used_a
        res[0]['positionA'] = posa
        res[0]['positionB'] = posa
        res[0]['N_spaxA'] = len(nsxA)
        res[0]['N_spaxB'] = len(nsxB)
        res[0]['meta'] = meta
        res[0]['object_spaxel_ids_A'] = nsxA
        res[0]['sky_spaxel_ids_A'] = kixa
        res[0]['object_spaxel_ids_B'] = nsxB
        res[0]['sky_spaxel_ids_B'] = kixb
        res[0]['sky_subtraction'] = False if stats['nosky'] else True
        res[0]['quality'] = quality
        # Calculate wavelength offsets
        coef = chebfit(np.arange(len(ll)), ll, 4)
        xs = np.arange(len(ll)+1)
        newll = chebval(xs, coef)
        # Store offsets
        res[0]['dlam'] = np.diff(newll)
        # Save the final spectrum
        np.save("sp_" + outname, res)
        print "Wrote sp_"+outname+".npy"


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=
        """Extract a spectrum from an image using a geometric cube solution.
Handles a single A image and A+B pair as well as flat extraction.
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--A', type=str, help='FITS A file')
    parser.add_argument('--B', type=str, help='FITS B file')
    parser.add_argument('fine', type=str, help='Numpy fine wavelength solution')
    parser.add_argument('--outname', type=str, help='Prefix output name')
    parser.add_argument('--std', type=str, help='Name of standard')
    parser.add_argument('--Aoffset', type=str,
                        help='Name of "A" flexure offset correction file')
    parser.add_argument('--radius_as', type=float,
                        help='Extraction radius in arcseconds', default=4)
    parser.add_argument('--flat_correction', type=str,
                        help='Name of flat field .npy file', default=None)
    parser.add_argument('--nosky', action="store_true", default=False,
                        help='No sky subtraction: only sum in aperture')
    parser.add_argument('--extflat', action="store_true", default=False,
                        help='Perform flat extraction')
    parser.add_argument('--specExtract', action="store_true", default=False,
                        help='Perform spectral extraction')
    parser.add_argument('--autoExtract', action="store_true", default=False,
                        help='Perform automatic extraction')
    parser.add_argument('--interact', action="store_true", default=False,
                        help='Interactively set the quality')

    args = parser.parse_args()

    print ""

    if args.outname is not None:
        args.outname = os.path.splitext(args.outname)[0]

    if args.flat_correction is not None:
        print "Using flat data in %s" % args.flat_correction
        flat = np.load(args.flat_correction)
    else:
        flat = None

    if args.A is not None and args.B is not None:
        print "A B Extraction to %s.npy" % args.outname
        handle_dual(args.A, args.B, args.fine, outname=args.outname,
                    offset=args.Aoffset, radius=args.radius_as,
                    flat_corrections=flat, nosky=args.nosky,
                    specExtract=args.specExtract,
                    interact=args.interact)

    elif args.A is not None:
        if args.std is None:
            if args.extflat:
                print "Flat Extraction to %s.npy" % args.outname
                handle_flat(args.A, args.fine, outname=args.outname)
            else:
                print "Single Extraction to %s.npy" % args.outname
                handle_single(args.A, args.fine, outname=args.outname,
                              offset=args.Aoffset, radius=args.radius_as,
                              flat_corrections=flat, nosky=args.nosky,
                              specExtract=args.specExtract,
                              autoExtract=args.autoExtract,
                              interact=args.interact)
        else:
            print "Standard Star Extraction to %s.npy" % args.outname
            star = Stds.Standards[args.std]
            handle_std(args.A, args.fine, outname=args.outname,
                       standard=star, offset=args.Aoffset,
                       flat_corrections=flat)

    else:
        print "I do not understand your intent, you must specify --A, at least"
