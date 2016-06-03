
import argparse
import numpy as np
import pylab as pl
import pyfits as pf

import NPK.Fit as NFit

from scipy.interpolate import interp1d

import SEDMr.Wavelength as Wavelength


def measure_flexure_x(cube, hdulist, drow=0., skylines=(557.0, 589.0),
                      lamstart=1000.0, lamratio=239./240., lamlen=250,
                      extract_width=3, skywidth=9, outfile='dX'):
    """Measures flexure in X direction, returns pixel offset

    Args:
        cube (extraction array): List of Extraction object, the fine loc +
            wave solution for each spectrum
        hdulist (pyfits obj): Pyfits object for the spectrum to measure
        drow (float): offset in rows for flexure (y-axis)
        skylines (float, float): The night skylines to centroid on in nm
        skywidth(float): Fit gaussian to the ROI of (skyline-skywidth to
            skyline+skywidth) in nm.
        extract_width(int): Number of pixels to extract spectrum around
        - See Wavelength.fiducial spectrum for following:
        lamstart (float): Wavelength to start the grid on, default 1000 nm
        lamratio (float): Resolution of sed machine
        lamlen (int): Length of spectrum
        outfile (string): output pdf plot showing fits to skyline(s)

    Returns:
        Offset number of pixels in X direction.
    """

    dat = hdulist[0].data

    spec_ixs = np.arange(500, 1200, 10)
    lamgrid = Wavelength.fiducial_spectrum(lamstart=lamstart,
                                           lamratio=lamratio, len=lamlen)

    specgrid = np.zeros((len(lamgrid), len(spec_ixs)))

    for i, ix in enumerate(spec_ixs):
        f = cube[ix]

        # bad fit
        if not f.ok:
            continue
        # noisy fit
        if f.lamnrms > 1:
            continue
        # short spectrum
        if f.xrange[1] - f.xrange[0] < 200:
            continue

        spec = np.zeros(f.xrange[1] - f.xrange[0])
        yfun = np.poly1d(f.poly)

        for jx, xpos in enumerate(np.arange(f.xrange[0], f.xrange[1])):
            ypos = yfun(xpos)

            try:
                spec[jx] = np.sum(dat[ypos-extract_width:ypos+extract_width,
                                  xpos])
            except:
                continue

        try:
            ll = f.get_lambda_nm()
        except:
            continue
        specfun = interp1d(ll, spec, bounds_error=False)
        specgrid[:, i] = specfun(lamgrid)

    skyspec = np.median(specgrid, axis=1)
    pl.step(lamgrid, skyspec, where='mid')
    pl.xlabel("Wavelength [nm]")
    pl.ylabel("Spec Irr [ph/10 m/nm]")

    sumoff = 0.
    sumscale = 0.
    for skyline in skylines:
        roi = (lamgrid > skyline-skywidth) & (lamgrid < skyline+skywidth)
        ffun = NFit.mpfit_residuals(NFit.gaussian4)
        parinfo = [
            {'value': np.max(skyspec[roi]), 'limited': [1, 0],
             'limits': [0, 0]},
            {'value': skyline},
            {'value': 3},
            {'value': np.min(skyspec[roi]), 'limited': [1, 0],
             'limits': [0, 0]}]
        fit = NFit.mpfit_do(ffun, lamgrid[roi], skyspec[roi], parinfo)
        if fit.status == 1:
            sumoff += (fit.params[1] - skyline) * fit.params[0]
            sumscale += fit.params[0]
        pl.plot(lamgrid[roi], NFit.gaussian4(fit.params, lamgrid[roi]))

        dxnm = fit.params[1] - skyline

        print "line: %6.1f (%6.1f), dX = %3.2f nm shift" % (skyline,
                                                            fit.params[0],
                                                            dxnm)

    if sumscale > 0.:
        dxnm = sumoff / sumscale
    else:
        print "Warning: no good skylines to fit!  Setting X offset to 0.0 nm"
        dxnm = 0.
    print "dX = %3.2f nm shift" % dxnm

    pl.title("dX = %3.2f nm shift, dY = %3.2f px shift" % (dxnm, drow))

    pl.savefig(outfile + ".pdf")

    return dxnm


def measure_flexure_y(cube, hdulist, profwidth=5, plot=False):

    dat = hdulist[0].data

    profs = []
    xx = np.arange(profwidth*2)

    for ix in np.arange(500, 1200, 10):
        f = cube[ix]
        profile = np.zeros(profwidth*2)

        # bad fit
        if not f.ok:
            continue
        # noisy fit
        if f.lamnrms > 1:
            continue
        # short spectrum
        if f.xrange[1] - f.xrange[0] < 200:
            continue

        yfun = np.poly1d(f.poly)
        for xpos in np.arange(f.xrange[0], f.xrange[1]):
            try:
                ypos = int(np.round(yfun(xpos)))
            except:
                continue
            try:
                profile += dat[ypos-profwidth:ypos+profwidth, xpos]
            except:
                continue

        profs.append(profile - np.min(profile))

    if plot:
        pl.figure(1)
    ffun = NFit.mpfit_residuals(NFit.gaussian4)
    parinfo = [{'value': 1}, {'value': profwidth}, {'value': 2},
               {'value': 0}, {'value': 0}]

    profposys = []
    for prof in profs:
        if plot:
            pl.step(xx, prof)
        parinfo[0]['value'] = np.max(prof)
        fit = NFit.mpfit_do(ffun, xx, prof, parinfo)
        if plot:
            pl.plot(xx, NFit.gaussian4(fit.params, xx))
        profposys.append(fit.params[1] - profwidth-1)
    if plot:
        pl.show()
    profposys = np.array(profposys)

    mn = np.mean(profposys)
    sd = np.std(profposys)
    ok = np.abs(profposys - mn)/sd < 3
    required_shift = np.mean(profposys[ok])
    print "dY = %3.2f pixel shift" % required_shift

    return required_shift


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

    args = parser.parse_args()

    fine, meta = np.load(args.fine)
    HDU = pf.open(args.infile)
    dy = measure_flexure_y(fine, HDU, profwidth=args.profwidth)
    dx = measure_flexure_x(fine, HDU, drow=dy,
                           skylines=args.skylines,
                           lamstart=args.lamstart,
                           lamratio=args.lamratio,
                           lamlen=args.lamlen,
                           extract_width=args.extract_width,
                           skywidth=args.skywidth,
                           outfile=args.outfile)

    res = [{'fine_name': args.fine,
            'infile_name': args.infile,
            'profwidth': args.profwidth,
            'skyline': args.skylines,
            'lamstart': args.lamstart,
            'lamratio': args.lamratio,
            'lamlen': args.lamlen,
            'extract_width': args.extract_width,
            'skywidth': args.skywidth,
            'dXnm': dx,
            'dYpix': dy}]

    np.save(args.outfile, res)
