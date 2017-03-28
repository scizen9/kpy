"""

 Finds SEDM spectra based on a sextractor segmentation map

    [c] 2014 Nick Konidaris

"""

import argparse
from collections import namedtuple
import numpy as np
import astropy.io.fits as pf
from multiprocessing import Pool
import sys
import warnings

from stsci.tools.gfit import gfit1d

import NPK.Bar as Bar

Coord = namedtuple("XY", "x y")


def spanrange(locs):
    """Returns a pair of coordinates where index values are true for a 2d array.

    Args:
        locs[d1,d2]: Boolean array of coordinates
    
    Returns:
        [[min x, min y], [max x, max y]]
        
    """

    global XX, YY, xi, yi

    # To avoid expensive creation
    try:
        XX
    except:
        xi = np.arange(locs.shape[0])
        yi = np.arange(locs.shape[1])
        XX, YY = np.meshgrid(xi, yi)

    mnx = np.min(XX[locs])
    mny = np.min(YY[locs])
    mxx = np.max(XX[locs])
    mxy = np.max(YY[locs])

    return [Coord(mnx, mny), Coord(mxx, mxy)]


def load_data(segmapfn, objfn):
    """Returns HDUs from files segmapfn and objfn.
    
    load_data is a thin wrapper around astropy.io.fits

    """

    try:
        segmap = pf.open(segmapfn)
    except Exception as e:
        print("Could not open %s: %s" % (segmapfn, e))

    try:
        obj = pf.open(objfn)
    except Exception as e:
        print("Could not open %s: %s" % (segmapfn, e))

    return segmap, obj


def find_segments_helper(seg_cnt):
    """Trace a single spaxel segment and return the solution

    Args:
        seg_cnt (int): the segment number

    Returns:
        dict: Returns a segmentation map Dictionary
        containing::

            {   "seg_cnt": Segment ID number,
                "xs": List of x positions of trace,
                "mean_ys": Measured Y position (average) of the trace,
                "coeff_ys": polyfit coefficients to the mean_ys,
                "trace_sigma": Width of trace in pixels (1 sigma),
                "ok": Trace has more than 50 pixels     }

    """
    # Global is for inter process communication
    global segdat, objdat, polyorder, n_done, update_rate

    n_done += 1
    if n_done % update_rate == 0:
        Bar.update()

    # Padding in pixels around trace in Y
    PAD = 2

    # Get this segment
    the_seg = (segdat == seg_cnt)

    # Test for tiny segment
    test = the_seg.nonzero()
    # trace data
    tr = {
        "seg_cnt": seg_cnt,
        "xs": np.array(np.nan),
        "mean_ys": np.array(np.nan),
        "coeff_ys": np.array(np.nan),
        "trace_sigma": np.array(np.nan),
        "ok": False,
        "bkg_ok": False
    }
    if len(test[0]) < 10 or len(test[1]) < 10:
        outstr = '\r%4.4i: %4.4i, TINY SEGMENT, %-5s' % (seg_cnt, len(test[0]),
                                                         tr['ok'])
        print(outstr)
        sys.stdout.flush()
    # We are OK
    else:
        # Get segment's geometry
        span = spanrange(the_seg)
        mnsr = np.max((span[0].y - PAD, 0))
        mxsr = np.min((span[1].y + PAD, segdat.shape[1] - 1))
        y_slc = slice(mnsr, mxsr)
        # y_off = (span[0].y + span[1].y) / 2.0

        # How long is it in X?
        n_el = span[1].x - span[0].x

        # How wide is it in Y?
        n_wd = span[1].y - span[0].y

        # Don't fit short traces, but flag them by setting "ok" to False
        if n_el < 50 or n_wd < 3:

            # Flag the sigma with zero
            sig = 0.
            tr["trace_sigma"] = sig

        # Trace is long enough and wide enough, so let's fit it!
        else:

            means = np.zeros(n_el)
            trace_profile = np.zeros(mxsr - mnsr)

            for i in range(n_el):
                with warnings.catch_warnings():
                    warnings.simplefilter("ignore", category=RuntimeWarning)
                    XX = i + span[0].x
                    profile = np.nanmedian(objdat[y_slc, XX - 3:XX + 3], 1)
                    profile -= np.min(profile)

                    trace_profile += profile

                    xs = np.arange(len(profile)) + span[0].y

                    means[i] = np.sum(xs * profile) / np.sum(profile) - PAD

            xs = np.arange(n_el) + span[0].x

            nans = ~np.isfinite(means)
            ok = np.isfinite(means)

            if np.count_nonzero(nans) < 5:
                poly = np.array(np.polyfit(xs[ok], means[ok], polyorder))
            else:
                poly = np.array(np.nan)

            trace_fit = gfit1d(trace_profile,
                               par=[0, len(trace_profile) / 2., 1.7], quiet=1)
            sig = np.abs(trace_fit.params[2])

            tr["xs"] = np.array(xs[ok])
            tr["mean_ys"] = np.array(means[ok])
            tr["coeff_ys"] = poly
            tr["trace_sigma"] = sig
            tr["bkg_ok"] = True
            tr["ok"] = (np.count_nonzero(nans) <= 0 & np.isfinite(poly[0]))

        # outstr = '\r%4.4i: %4.4i, fwhm=%3.2f pix, %-5s' % (seg_cnt, n_el,
        #                                                    sig * 2.355,
        #                                                    tr['ok'])
        # print(outstr, end="")
        # sys.stdout.flush()

    return tr


def find_segments(segmap=None, obj=None, order=2):
    """Find the segments in obj by tracing ridgelines identified in 
    the segmentation map.

    Args:
        segmap ([int,int]): Segmentation map image, first segment is
                1 .. max(segmap)
        obj ([float,float]): Image data to trace over
        order (int): The order of the polynomial used in coeff_ys

    Returns:
        list of dict: Returns a list with length equal to
        the max(segmap) segmentation map Dictionaries containing::

            {   "seg_cnt": Segment ID number,
                "xs": List of x positions of trace,
                "mean_ys": Measured Y position (average) of the trace,
                "coeff_ys": polyfit coefficients to the mean_ys,
                "trace_sigma": Width of trace in pixels (1 sigma),
                "ok": Trace has more than 50 pixels     }
        
    """
    global segdat, objdat, polyorder, n_done, update_rate

    segdat = segmap[0].data
    objdat = obj[0].data
    polyorder = order

    # First segment begins at 1
    segrange = range(1, max(segdat.flatten())+1)

    n_done = 0
    update_rate = int(len(segrange) / Bar.setup(toolbar_width=74)) + 1
    p = Pool(8)
    traces = p.map(find_segments_helper, segrange)
    p.close()
    Bar.done(mapped=True)
    # print("")

    return traces


def ds9_line_str(startx, starty, endx, endy):
    """Returns a string representing a line in ds9"""

    return "line(%s,%s,%s,%s) # line=0 0\n" % (startx, starty, endx, endy)


def write_reports(segments, outname):
    ds9 = "global color=blue\nphysical\n"
    for seg in segments:
        if not seg['ok']:
            continue

        ff = np.poly1d(seg['coeff_ys'])
        for dx in [0, 100, 200]:
            X1 = seg['xs'][0] + dx
            Y1 = ff(X1) + 1
            X2 = X1 + 100
            Y2 = ff(X2) + 1

            if np.isfinite(X1) and np.isfinite(X2) and np.isfinite(Y1) and \
                    np.isfinite(Y2):
                ds9 += ds9_line_str(X1, Y1, X2, Y2)

    fname = "ds9_%s.reg" % outname
    f = open(fname, 'w')
    f.write(ds9)
    f.close()


def write_segments(segments):
    np.save(objfn + "_segments", segments)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=
                                     """
    FindSpectra.py
    """)

    parser.add_argument("segmap_fits", type=str,
                        help="Segmentation map .fits file name")
    parser.add_argument("objfn", type=str, help="Dome flat")
    parser.add_argument("outname", type=str, help="Output file name prefix")
    parser.add_argument("--order", type=int, help="Polynomial order",
                        default=2)

    args = parser.parse_args()
    segmapfn = args.segmap_fits
    objfn = args.objfn
    outname = args.outname
    order = args.order

    segmap, obj = load_data(segmapfn, objfn)
    segments = find_segments(segmap, obj, order=order)
    write_segments(segments)
    write_reports(segments, outname)
