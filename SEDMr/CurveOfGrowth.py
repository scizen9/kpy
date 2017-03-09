
import argparse
import os
import numpy as np
import pylab as pl
import itertools

from astropy.coordinates import Angle
from numpy.polynomial.chebyshev import chebval
from scipy.interpolate import interp1d
import SEDMr.Wavelength as Wavelength
import SEDMr.Spectra as SedSpec
import SEDMr.Version as Version

# Nadia imports
from scipy.interpolate import griddata
import scipy.optimize as opt

import skimage.feature as feature

drp_ver = Version.ifu_drp_version()


def plot_drp_ver():
    ax = pl.gca()
    ax.annotate('DRP: ' + drp_ver, xy=(0.0, 0.01), xytext=(0, 0),
                xycoords=('axes fraction', 'figure fraction'),
                textcoords='offset points', size=6,
                ha='center', va='bottom')


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
                               airmass=1.0, sigfac=3.0, plotobj=False):
    """ 
    Returns index of spectra picked by Guassian fit.
    
    NOTE: Index is counted against the array, not seg_id
    """

    status = 0

    pl.ioff()

    kt = SedSpec.Spectra(spectra)

    # Get X,Y positions (arcsec) and summed values between lmin and lmax
    xs, ys, vs = kt.to_xyv(lmin=lmin, lmax=lmax)

    xi = np.linspace(np.nanmin(xs), np.nanmax(xs), 200)
    yi = np.linspace(np.nanmin(ys), np.nanmax(ys), 200)

    x, y = np.mgrid[np.nanmin(xs):np.nanmax(xs):200j,
                    np.nanmin(ys):np.nanmax(ys):200j]

    points = zip(xs, ys)
    values = vs

    # Create image, print(stats)
    grid_vs = griddata(points, values, (x, y), method='linear')
    grid_vs[np.isnan(grid_vs)] = np.nanmean(grid_vs)
    grid_med = np.nanmedian(grid_vs)
    print("grid_vs min, max, mean, median: %f, %f, %f, %f\n" %
          (np.nanmin(grid_vs), np.nanmax(grid_vs),
           np.nanmean(grid_vs), grid_med))

    # Find features in image
    blobs = feature.blob_dog(grid_vs-grid_med, min_sigma=10, max_sigma=20,
                             threshold=100.0)
    print("Found %d blobs" % len(blobs))

    goodblob = 0

    # Loop over found blobs
    objs = []
    for blob in blobs:
        # Extract blob properties
        bx, by, br = blob

        bx = int(bx)
        by = int(by)

        # How bright is this blob?
        gv = grid_vs[bx, by]-grid_med
        # Exclude edge blobs and faint blobs
        if 0 < bx < 199 and 0 < by < 199 and gv > 100.:
            goodblob += 1
            print("%3d, z, x, y, dra, ddec: %8.1f, %5d, %5d, %6.2f, %6.2f" %
                  (goodblob, gv, bx, by, xi[bx], yi[by]))
            objs.append((gv, xi[bx], yi[by], goodblob))

    print("Found %d good objects" % len(objs))
    if len(objs) <= 0:
        objs = [(1000., 0., 0., goodblob)]
    # Make sure the brightest object is last
    objs.sort()

    # Perform 2-D Gaussian fit of good (real) objects
    sigma_x = 2.0
    sigma_y = 2.0
    for obj in objs:
        # Fill initial fit params
        amplitude = obj[0]
        xo = obj[1]
        yo = obj[2]
        objno = obj[3]

        print("\nFitting object %d" % objno)
        print("initial guess : z,a,b,x,y,theta:"
              " %9.1f, %6.2f, %6.2f, %6.2f, %6.2f, %7.2f" %
              (amplitude, sigma_x, sigma_y, xo, yo, 0.))

        # create initial data
        initial_guess = (amplitude, xo, yo, sigma_x, sigma_y, 0, grid_med)

        try:
            popt, pcov = opt.curve_fit(gaussian_2d, (x, y),
                                       grid_vs.flatten(), p0=initial_guess)
        except RuntimeError:
            print("ERROR: unable to fit Gaussian")
            print("Using initial guess")
            status = 3
            popt = initial_guess
        # Fitted position
        xc = popt[1]
        yc = popt[2]
        if xc < -30. or xc > 30. or yc < -30. or yc > 30.:
            print("ERROR: X,Y out of bounds: %f, %f" % (xc, yc))
            print("Using initial position: %f, %f" % (xo, yo))
            xc = xo
            yc = yo
            status = 1
        pos = (xc, yc)
        # Fitted 3-sigma extent
        a = popt[3]*sigfac
        b = popt[4]*sigfac
        if a > 14. or b > 14. or a <= 0. or b <= 0.:
            print("ERROR: A,B out of bounds: %f, %f" % (a, b))
            print("Using initial values: %f, %f" % (sigma_x, sigma_y))
            a = sigma_x
            b = sigma_y
            status = 2
        theta = popt[5]
        z = popt[0]

        # report position and shape
        ellipse = (a, b, xc, yc, theta * (180. / np.pi))
        print("PSF FIT on IFU: z,a,b,x,y,theta:"
              " %9.1f, %6.2f, %6.2f, %6.2f, %6.2f, %7.2f\n" %
              (z, a, b, xc, yc, theta*180./np.pi))

        positions = [pos]

        # Gather spaxels
        all_kix = []
        for the_pos in positions:
            all_kix.append(list(find_positions_ellipse(kt.KT.data, xc, yc, a, b,
                                                       -theta)))

        all_kix = list(itertools.chain(*all_kix))
        kix = list(set(all_kix))
        print("found this many spaxels: %d" % len(kix))

    if status == 0 and goodblob == 0:
        print("ERROR: no good objects found in image")
        status = 4

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
        print("Number of starting pure sky spaxels is %d" % len(skys))
    else:
        print("ERROR: no sky spaxels in this image")

    newspec = [spectra[i] for i in skys]
    kt = SedSpec.Spectra(newspec)

    xs, ys, vs = kt.to_xyv(lmin=lmin, lmax=lmax)
    vmdn = np.nanmedian(vs)

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

        if np.nanmedian(el.spec[ok]) > hi_thresh:
            skys.remove(s)
            n_hi_rem += 1

        if np.nanmedian(el.spec[ok]) < lo_thresh:
            skys.remove(s)
            n_lo_rem += 1

        n_tot += 1

    n_tot -= n_hi_rem + n_lo_rem
    print("Removed %d high sky spaxels and %d low sky spaxels leaving %d "
          "remaining spaxels" % (n_hi_rem, n_lo_rem, n_tot))

    return skys


def c_to_nm(coefficients, pix, offset=0.):

    t = coefficients[:]
    t[0] += offset
    return chebval(pix, t)


def interp_spectra(all_spectra, six, onto=None, sky=False):
    """Interp spectra onto common grid

    Args:
        all_spectra:
        six:
        onto:
        sky:
    """

    l_grid = onto
    s_grid = []
    lamcoeff = None
    # for ix,spectrum in enumerate(all_spectra):
    for ix in six:
        spectrum = all_spectra[ix]

        if sky and all_spectra[ix].is_obj:
            continue

        l, s = spectrum.get_counts(the_spec='specw')
        pix = np.arange(*spectrum.xrange)

        # check for saturated traces
        if np.max(s) > 1000000:
            print("saturated extraction: %d with max of %d, skipping" %
                  (ix, np.max(s)))
            continue

        # This is correct: preference to lamcoeff over mdn_coeff
        if spectrum.lamcoeff is not None:
            cs = spectrum.lamcoeff
        else:
            cs = spectrum.mdn_coeff

        # get wavelengths for spectrum
        l = c_to_nm(cs, pix)

        # skip short spectra (on or near edge of IFU)
        if l.max() - l.min() < 300:
            continue

        # Check if our wavelength grid is defined,
        if l_grid is None:
            # use the first set of wavelengths and store
            l_grid = l
            fl = s
            lamcoeff = spectrum.lamcoeff
        else:
            # Interpolate onto our wavelength grid and store
            fun = interp1d(l, s, bounds_error=False, fill_value=0)
            fl = fun(l_grid)

        s_grid.append(fl)

    # average of all spectra selected
    medspec = np.nanmean(s_grid, axis=0)

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

    return result


def make_cog(outname, lmin=650., lmax=700., sigfac=7., interact=False):
    """Loads IFU frame "imfile" and extracts spectra using "fine".

    Args:
        outname (string): filename to write results to
        lmin (float): lower wavelength limit for image generation
        lmax (float): upper wavelength limit for image generation
        sigfac (float): sigma multiplier for Gaussian extent of aperture
        interact (bool): set for interactive plotting

    Returns:
        None

    Raises:
        None

    Note:
        The extracted spectrum dictionary is written to ``outname``.npy::

            'doc': docstring description
            'ph_10m_nm': Flux in photon/10min/nm integrated
            'spectra': Individual spaxel spectra
            'coefficients': Coefficients of wavelength fit
            'nm': Wavelengths in nm


    """

    # The spaxel extraction file must already exist, so load extractions in
    if os.path.isfile(outname+".npy"):
        print("USING extractions in %s.npy!" % outname)
        ex, meta = np.load(outname+".npy")
    # No file found
    else:
        print("File not found: %s.npy" % outname)
        return

    # Get the object name of record
    objname = meta['header']['OBJECT'].split()[0]

    # Automatic extraction using Gaussian fit for Standard Stars
    sixa, posa, adcpos, ellipse, status = \
        identify_spectra_gauss_fit(ex,
                                   prlltc=Angle(meta['PRLLTC'], unit='deg'),
                                   lmin=lmin, lmax=lmax,
                                   airmass=meta['airmass'],
                                   sigfac=sigfac)
    radius_used = ellipse[0] * 0.5

    if status > 0:
        quality = 4  # something went wrong
    else:
        quality = 1  # by definition

    # Stats for automatic fit
    stats = {"nosky": False, "scaled": False,
             "lmin": lmin, "lmax": lmax,
             "cmin": None, "cmax": None}

    # Use all sky spaxels in image
    kixa = identify_sky_spectra(ex, posa, ellipse=ellipse)

    for ix in sixa:
        ex[ix].is_obj = True
    for ix in kixa:
        ex[ix].is_sky = True

    # Get sky spectrum
    skya = interp_spectra(ex, kixa, sky=True)

    # Define our standard wavelength grid
    ll = Wavelength.fiducial_spectrum()
    # Resample sky onto standard wavelength grid
    sky_a = interp1d(skya[0]['nm'], skya[0]['ph_10m_nm'],
                     bounds_error=False)

    # Set up curve of growth
    kt = SedSpec.Spectra(ex)
    elrat = ellipse[1] / ellipse[0]
    xc = ellipse[2]
    yc = ellipse[3]
    theta = ellipse[4]*np.pi/180.

    # Set up plot
    pl.figure(1)
    pl.clf()
    pl.title(objname)

    xs = range(20)
    rs = list(np.linspace(0, ellipse[0], 20))
    rs.reverse()
    print("loop down to 1 from %.2f asec" % ellipse[0])

    c1 = []
    c2 = []
    c3 = []
    c4 = []
    c5 = []
    px = []
    resout = []

    # Loop over semi-major axes
    for ix in xs:
        if rs[ix] <= 0.:
            continue
        a = rs[ix]
        b = rs[ix] * elrat
        kix = find_positions_ellipse(kt.KT.data, xc, yc, a, b, -theta)
        sixa = kt.good_positions[kix]
        print("%02d found %04d spaxels with a = %.3f" % (ix, len(sixa), a))

        if len(sixa) > 0:
            # get the summed spectrum over the selected spaxels
            resa = interp_spectra(ex, sixa)

            # Copy and resample object spectrum onto standard wavelength grid
            res = {
                "doc": resa[0]["doc"],
                "ph_10m_nm": np.copy(resa[0]["ph_10m_nm"]),
                "spectra": np.copy(resa[0]["spectra"]),
                "coefficients": np.copy(resa[0]["coefficients"]),
                "nm": np.copy(resa[0]["ph_10m_nm"]), 'nm': np.copy(ll)
                }
            fl = interp1d(resa[0]['nm'], resa[0]['ph_10m_nm'],
                          bounds_error=False)
            # Calculate output corrected spectrum
            # Account for sky and aperture
            res['ph_10m_nm'] = (fl(ll)-sky_a(ll)) * len(sixa)
            cog = res['ph_10m_nm']
            # 400 - 500 nm
            f1 = np.nanmean(cog[(ll > 400) * (ll < 500)])
            c1.append(f1)
            # 500 - 600 nm
            f2 = np.nanmean(cog[(ll > 500) * (ll < 600)])
            c2.append(f2)
            # 600 - 700 nm
            f3 = np.nanmean(cog[(ll > 600) * (ll < 700)])
            c3.append(f3)
            # 700 - 800 nm
            f4 = np.nanmean(cog[(ll > 700) * (ll < 800)])
            c4.append(f4)
            # 800 - 900 nm
            f5 = np.nanmean(cog[(ll > 800) * (ll < 900)])
            c5.append(f5)
            # Semi-major axis in arcsec
            px.append(a)
            # Output results
            res['a'] = a
            res['ix'] = ix
            res['Nspax'] = len(sixa)
            res['fl_400_500nm'] = f1
            res['fl_500_600nm'] = f2
            res['fl_600_700nm'] = f3
            res['fl_700_800nm'] = f4
            res['fl_800_900nm'] = f5
            # Store this iteration
            resout.append(res)

    maxph = np.nanmax([c1, c2, c3, c4, c5])

    hmx = np.array(px)

    half = (c3 > maxph/2.)
    hmx = hmx[half]
    hmx = hmx[-1]

    c1 /= maxph
    c2 /= maxph
    c3 /= maxph
    c4 /= maxph
    c5 /= maxph

    pl.plot(px, c1, label='400-500 nm')
    pl.plot(px, c2, label='500-600 nm')
    pl.plot(px, c3, label='600-700 nm')
    pl.plot(px, c4, label='700-800 nm')
    pl.plot(px, c5, label='800-900 nm')
    pl.plot([hmx, hmx], [-0.05, 1.05], ls='--', c='black', label='HalfMax')
    pl.ylim(-0.05, 1.05)
    pl.xlabel('Semi-major axis (arcsec)')
    pl.ylabel('Relative Irr')
    pl.legend()
    plot_drp_ver()
    if interact:
        pl.show()
    else:
        pl.savefig('cog_' + outname + '.pdf')
        print("Wrote cog_"+outname+".pdf")

    np.save("cog_" + outname, resout)
    print("Wrote cog_"+outname+".npy")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description=
        """Plot a curve of growth for the brightest object in the extractions.
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('--inname', type=str, help='Input extractions name')
    parser.add_argument('--sigfac', type=float,
                        help='Gaussian Sigma multiplier for maximum extent of'
                             ' aperture', default=7.0)
    parser.add_argument('--interact', action='store_true', default=False,
                        help='Interactive plotting')

    args = parser.parse_args()

    print("")

    make_cog(args.inname, sigfac=args.sigfac, interact=args.interact)
