
import argparse
import numpy as np
import pylab as pl
import astropy.io.fits as pf

import NPK.Fit as NFit

from scipy.interpolate import interp1d

import SEDMr.Wavelength as Wavelength


def measure_flexure_x(cube, hdulist, drow=0., skylines=(557.0, 589.0),
                      lamstart=1000.0, lamratio=239./240., lamlen=250,
                      extract_width=3, skywidth=9, outfile='dX', plot=False):
    """Measures flexure in X direction, returns pixel offset

    Args:
        cube (extraction array): List of Extraction object, the fine loc +
            wave solution for each spectrum
        hdulist (astropy.io.fits obj): Pyfits object for the spectrum to measure
        drow (float): offset in rows for flexure (y-axis)
        skylines (float, float): The night skylines to centroid on in nm

        - See Wavelength.fiducial spectrum for following:
        lamstart (float): Wavelength to start the grid on, default 1000 nm
        lamratio (float): Resolution of sed machine
        lamlen (int): Length of spectrum

        extract_width(int): Number of pixels to extract spectrum around

        skywidth(float): Fit gaussian to the ROI of skyline+-skywidth in nm.
        outfile (string): output pdf plot showing fits to skyline(s)
        plot (bool): set to True to show plot, else plot to pdf file

    Returns:
        Offset number of pixels in X direction.
    """

    # read in image
    dat = hdulist[0].data
    # select 70 representative spaxels
    spec_ixs = np.arange(500, 1200, 10)
    # fiducial wavelength grid
    lamgrid = Wavelength.fiducial_spectrum(lamstart=lamstart,
                                           lamratio=lamratio, npx=lamlen)
    # initialize grid of spaxel spectra
    specgrid = np.zeros((len(lamgrid), len(spec_ixs)))

    # loop over selected spaxels
    for i, ix in enumerate(spec_ixs):
        # get spaxel
        f = cube[ix]

        # skip baddies
        # bad fit
        if not f.ok:
            continue
        # noisy fit
        if f.lamnrms > 1:
            continue
        # short spectrum
        if f.xrange[1] - f.xrange[0] < 200:
            continue

        # set up spectral vector for this spaxel
        spec = np.zeros(f.xrange[1] - f.xrange[0])
        yfun = np.poly1d(f.poly)

        # loop over x positions in spaxel
        for jx, xpos in enumerate(np.arange(f.xrange[0], f.xrange[1])):

            # get y position on image
            ypos = int(yfun(xpos))

            # extract spectrum from image
            try:
                spec[jx] = np.sum(dat[ypos-extract_width:ypos+extract_width,
                                  xpos])
            except:
                print "Warning: no sum for sky spectrum %d at %d" % (jx, xpos)
                continue

        # get wavelengths of spaxel
        try:
            ll = f.get_lambda_nm()
        except:
            continue

        # resample spectrum on fiducial wavelength grid
        specfun = interp1d(ll, spec, bounds_error=False)
        # insert into grid of spaxel spectra
        specgrid[:, i] = specfun(lamgrid)

    # create a median spectrum from spaxel grid
    # taking a median minimizes impact of objects in sample
    skyspec = np.median(specgrid, axis=1)

    # plot resulting sky spectrum
    pl.step(lamgrid, skyspec, where='mid')
    pl.xlabel("Wavelength [nm]")
    pl.ylabel("Spec Irr [ph/10 m/nm]")
    legend = ["Sky"]

    # accumulate average offsets from known sky lines
    sumoff = 0.
    sumscale = 0.
    minsig = 10000.

    # loop over input sky lines
    for skyline in skylines:

        # extract a wavelength window around sky line
        roi = (lamgrid > skyline-skywidth) & (lamgrid < skyline+skywidth)
        # prepare for Gaussian fit
        ffun = NFit.mpfit_residuals(NFit.gaussian4)
        # initial setup of fit
        parinfo = [
            {'value': np.max(skyspec[roi]), 'limited': [1, 0],
             'limits': [0, 0]},
            {'value': skyline},
            {'value': 3},
            {'value': np.min(skyspec[roi]), 'limited': [1, 0],
             'limits': [0, 0]}]
        # do the fit
        fit = NFit.mpfit_do(ffun, lamgrid[roi], skyspec[roi], parinfo)
        # did the fit succeed?
        if fit.status == 1 and fit.params[2] > 0.:
            off = fit.params[1] - skyline
            sumoff += off * fit.params[0]
            sumscale += fit.params[0]
            if fit.params[2] < minsig:
                minsig = fit.params[2]
            pl.plot(lamgrid[roi], NFit.gaussian4(fit.params, lamgrid[roi]))
            dxnm = fit.params[1] - skyline
            legend.append("%.1f, %.2f" % (skyline, off))
        else:
            dxnm = 0.

        print("line = %6.1f (%6.1f), FWHM = %.2f nm, status = %d, dX = %3.2f nm shift" %
              (skyline, fit.params[0], fit.params[2]*2.354, fit.status, dxnm))

    if sumscale > 0.:
        dxnm = sumoff / sumscale
    else:
        print "Warning: no good skylines to fit!  Setting X offset to 0.0 nm"
        dxnm = 0.
        minsig = 0.
    print "dX = %3.2f nm shift" % dxnm

    pl.title("dX = %3.2f nm shift, dY = %3.2f px shift" % (dxnm, drow))
    pl.legend(legend)

    if plot:
        pl.show()
    else:
        pl.savefig(outfile + ".pdf")

    # return offset and best FWHM
    return dxnm, minsig*2.354


def measure_flexure_y(cube, hdulist, profwidth=5, plot=False):

    # get image to measure
    dat = hdulist[0].data

    # store profiles
    profs = []

    # x values of profiles
    xx = np.arange(profwidth*2)

    # get a sample (70) of spaxels
    for ix in np.arange(500, 1200, 10):

        # get specific spaxel
        f = cube[ix]

        # set up profile vector for this spaxel
        profile = np.zeros(profwidth*2)

        # weed out baddies
        # bad fit
        if not f.ok:
            continue
        # noisy fit
        if f.lamnrms > 1:
            continue
        # short spectrum
        if f.xrange[1] - f.xrange[0] < 200:
            continue

        # get trace of spaxel on image in y
        yfun = np.poly1d(f.poly)

        # loop over x positions of this spaxel
        for xpos in np.arange(f.xrange[0], f.xrange[1]):

            # get the y position at the given xpos
            try:
                ypos = int(np.round(yfun(xpos)))
            except:
                continue

            # sum profile over all xpos
            try:
                profile += dat[ypos-profwidth:ypos+profwidth, xpos]
            except:
                continue

        # append profile after subtracting minimum
        profs.append(profile - np.min(profile))

    if plot:
        pl.figure(1)

    # set up Gaussian fit to profiles
    ffun = NFit.mpfit_residuals(NFit.gaussian4)
    # initial guess
    parinfo = [{'value': 1}, {'value': profwidth}, {'value': 2},
               {'value': 0}, {'value': 0}]

    # y positions of profiles
    profposys = []
    profwidys = []
    for prof in profs:
        # plot input profile
        if plot:
            pl.plot(xx, prof, 'ro')
        # update initial guess
        parinfo[0]['value'] = np.max(prof)
        # fit Gaussian
        fit = NFit.mpfit_do(ffun, xx, prof, parinfo)
        # overplot fit
        if plot:
            pl.plot(xx, NFit.gaussian4(fit.params, xx))
        # x offset between nominal trace position (profwidth-1) and
        # Gaussian fit position (fit.params[1])
        profposys.append(fit.params[1] - profwidth-1)
        profwidys.append(fit.params[2])

    profposys = np.array(profposys)
    profwidys = np.array(profwidys)

    # get statistics
    mn = np.mean(profposys)
    sd = np.std(profposys)
    # clean 3 sigma outliers
    ok = np.abs(profposys - mn)/sd < 3
    # final shift
    required_shift = np.mean(profposys[ok])
    # now do width
    mn = np.mean(profwidys)
    sd = np.std(profwidys)
    # clean 3 sigma outliers
    ok = np.abs(profwidys - mn)/sd < 3
    average_width = np.mean(profwidys[ok]) * 2.354
    print "yFWHM = %5.2f pixels" % average_width
    print "dY = %3.2f pixel shift" % required_shift

    if plot:
        px = (profwidth+1) + required_shift
        pl.plot([px, px], [0, np.max(profs)], '--')
        pl.show()

    return required_shift, average_width


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
        """Measure the flexure in x [nm] and y [px].

        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('fine', type=str, help='Fine correction path')
    parser.add_argument('infile', type=str, help='Path to FITS file to refit')
    parser.add_argument('--profwidth', type=float,
                        help='Profile width to extract for Y flexure',
                        default=5)
    parser.add_argument('--skylines', type=float,
                        help='skyline positions in nm to measure X flexure',
                        default=(557.7,589.0))
    parser.add_argument('--lamstart', type=float,
                        help='Wavelength to start interpolating grid ',
                        default=1000.0)
    parser.add_argument('--lamratio', type=float,
                        help='Wavelength resolution for interpolating grid',
                        default=239.0/240.0)
    parser.add_argument('--lamlen', type=int,
                        help='Wavelength grid length', default=250)
    parser.add_argument('--extract_width', type=int,
                        help='Extraction width in Y (for X flexure)', default=3)
    parser.add_argument('--skywidth', type=float,
                        help='Wavelength search window', default=9)
    parser.add_argument('--outfile', type=str, help='Output filename',
                        default="flexure.npy")
    parser.add_argument('--plot', action="store_true", default=False,
                        help='Plot results')

    args = parser.parse_args()

    fine, meta = np.load(args.fine)
    HDU = pf.open(args.infile)
    dy, ywid = measure_flexure_y(fine, HDU, profwidth=args.profwidth,
                                 plot=args.plot)
    dx, xwid = measure_flexure_x(fine, HDU, drow=dy,
                                 skylines=args.skylines,
                                 lamstart=args.lamstart,
                                 lamratio=args.lamratio,
                                 lamlen=args.lamlen,
                                 extract_width=args.extract_width,
                                 skywidth=args.skywidth,
                                 outfile=args.outfile,
                                 plot=args.plot)

    res = [{'fine_name': args.fine,
            'infile_name': args.infile,
            'profwidth': args.profwidth,
            'skyline': args.skylines,
            'lamstart': args.lamstart,
            'lamratio': args.lamratio,
            'lamlen': args.lamlen,
            'extract_width': args.extract_width,
            'skywidth': args.skywidth,
            'yfwhm': ywid,
            'xfwhm': xwid,
            'dXnm': dx,
            'dYpix': dy}]

    np.save(args.outfile, res)
