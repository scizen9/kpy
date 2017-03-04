"""Convert an extracted file into a data cube.

    STEP is either:
        * make: To create the data cube (once per night)
        * extract: To extract the cube (one for each observation)

Functions:
    * :func:`QR_to_img` convert cube to an image
    * :func:`extraction_to_cube` calculates sky positions
    * :func:`reject_outliers` trims data values of outliers

References:

    See figures here: http://www.astro.caltech.edu/sedm/_images/HexCoords.png

    Axial coordinates of a cube::

            (q,r-1)   (q+1,r-1)
        (q-1, r)   (q, r)   (q+1, r)
            (q-1,r+1)   (q, r+1)

    Even Q coordinates look like::
    
            0,0     1,0     2,0     3,0
        0,1     1,1     2,1     3,1     4,1
            0,2     1,2     2,2     3,2

Note:
    This is used as a python script as follows::

        usage: Cube.py [-h] [--step STEP] [--outname OUTNAME] extracted

        positional arguments:
          extracted          Extracted file

        optional arguments:
          -h, --help         show this help message and exit
          --step STEP        [make|extract|dump]
          --outname OUTNAME  Output cube name

"""

import os
import argparse
import numpy as np
import astropy.io.fits as pf
from scipy.spatial import KDTree 

from numpy.polynomial.chebyshev import chebval
from scipy.interpolate import interp1d
import SEDMr.Wavelength as Wavelength

import sys

sys.setrecursionlimit(10000)

ncol = 0
scale = 1.0
H2P = np.array([[np.sqrt(3), np.sqrt(3)/2], [0, 3/2.]]) * scale
P2H = np.array([[np.sqrt(3)/3, -1/3.], [0, 2/3.]]) / scale

# Rotation matrix
theta = np.deg2rad(-37+13.5)
ROT = np.array([[np.cos(theta), -np.sin(theta)], 
                [np.sin(theta),  np.cos(theta)]])


def qr_to_img(exts, size=4, outname="cube.fits"):
    """Convert a data cube to a fits image

    Args:
        exts (list of Extraction): extractions to convert (see Extraction.py)
        size (int): expansion factor, defaults to 4
        outname (str): output fits file name, defaults to cube.fits

    Returns:
        None

    """
    
    xs = np.array([ex.X_as for ex in exts], dtype=np.float)
    ys = np.array([ex.Y_as for ex in exts], dtype=np.float)

    minx = size * np.nanmin(xs)
    miny = size * np.nanmin(ys)
    maxx = size * np.nanmax(xs)
    maxy = size * np.nanmax(ys)

    dx = (maxx-minx)/.25
    dy = (maxy-miny)/.25
    l_grid = Wavelength.fiducial_spectrum()
    l_grid = l_grid[::-1]
    dl_grid = np.diff(l_grid)
    l_grid = l_grid[1:]

    img = np.zeros((dx, dy, len(l_grid)/2))
    img[:] = np.nan

    xsz = img.shape[0]/2
    ysz = img.shape[1]/2

    allspec = np.zeros((len(exts), len(l_grid)/2))
    for cnt, ex in enumerate(exts):
        if ex.xrange is None:
            continue
        if ex.exptime is None:
            ex.exptime = 1
        if ex.lamcoeff is None:
            continue

        ix = np.arange(*ex.xrange)
        l = chebval(ix, ex.lamcoeff)
        s = ex.specw

        f = interp1d(l, s, fill_value=np.nan, bounds_error=False)
        fi = f(l_grid) / dl_grid
        fi = fi[0:len(fi):2] + fi[1:len(fi):2]

        allspec[cnt, :] = fi

        x = (ex.X_as - minx)/0.25
        y = (ex.Y_as - miny)/0.25

        try: 
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    img[x+dx, y+dy, :] = fi
        except:
            pass

        # outstr = "\rX = %+10.5f, Y = %+10.5f" % (x[0], y[0])
        # print(outstr, end="")
        sys.stdout.flush()

    back = np.nanmedian(allspec, 0)

    if 'fits' not in outname:
        outname += '.fits'

    ff = pf.PrimaryHDU(img.T)
    ff.writeto(outname)

    for cnt, ex in enumerate(exts):
        if ex.xrange is None:
            continue
        if ex.exptime is None:
            ex.exptime = 1
        if ex.lamcoeff is None:
            continue

        ix = np.arange(*ex.xrange)
        l = chebval(ix, ex.lamcoeff)
        s = ex.specw

        f = interp1d(l, s, fill_value=np.nan, bounds_error=False)
        fi = f(l_grid)/dl_grid 
        fi = fi[0:len(fi):2] + fi[1:len(fi):2] - back

        try:
            x = np.round(size * np.sqrt(3.) * (ex.Q_ix + ex.R_ix / 2)) + xsz
            y = np.round(size * 3. / 2. * ex.R_ix) + ysz
        except:
            continue
        try: 
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    img[x+dx, y+dy, :] = fi
        except:
            pass

    ff = pf.PrimaryHDU(img.T)
    ff.writeto("bs_" + outname)


def extraction_to_cube(exts, outname="G.npy"):
    """ Convert the extraction to sky coordinates

    Args:
        exts (list of Extraction): The list of extractions (see Extraction.py)
        outname (str): The file created with the new extraction

    Returns:
        None

    Note:
        The new data cube has the following coordinate positions populated:

        * X_as: X position in arcsec
        * Y_as: Y position in arcsec
        * Z_as: Z position in arcsec
        (the Z coordinate runs 45 degrees to X and is not a `3rd` dimension).
        * Q_ix: The axial Q coordinate in integral units
        * R_ix: The axial R coordinate in integral units

        The relationship of Q/R to X/Y is defined through the pixel mapping
        matrix times the rotation matrix:

        * P2H = np.array([[np.sqrt(3)/3, -1/3.], [0, 2/3.]])
        * rot (22 degree)

    """

    global ncol

    # start with our x,y pixel positions
    xs = [None] * len(exts)
    ys = [None] * len(exts)

    # these are the hex axial coords
    for ex in exts:
        ex.Q_ix = None
        ex.R_ix = None
    
    # counters for which wavelength solution is used:
    # mdn = median, lam = lambda, tot = total used
    n_mdn = 0
    n_lam = 0
    n_tot = 0

    for ix, ex in enumerate(exts):
        # The X/Y location of a lenslet is based on its
        # trace Y position and where the fiducial wavelength
        # is expected in the X direction.

        # skip bad solutions
        if not ex.ok:
            continue

        # initialize
        xs[ix] = -999
        ys[ix] = -999

        # lamcoeff has precedence
        if ex.lamcoeff is not None:
            coeff = ex.lamcoeff
            n_lam += 1
        elif ex.mdn_coeff is not None:
            coeff = ex.mdn_coeff
            n_mdn += 1
        else:
            continue

        n_tot += 1

        # get pixel and wavelength vectors
        ixs = np.arange(*ex.xrange)
        ll = chebval(ixs, coeff)
        
        # fill in the x,y pix positions
        xs[ix] = np.nanmin(ex.xrange) + ex.xrefpix
        ys[ix] = np.nanmean(ex.yrange)
        ex.X_pix = xs[ix]
        ex.Y_pix = ys[ix]
        ex.Xccd_as = (1024.0 - xs[ix]) * 0.0128
        ex.Yccd_as = (1024.0 - ys[ix]) * 0.0128
        if ex.xrefpix is not None:
            if ex.xrefpix < len(ll):
                ex.xreflam = ll[ex.xrefpix]
    # End loop over all extractions

    # Record wavelength mean and report stats
    xreflams = reject_outliers(
            np.array([ex.xreflam for ex in exts], dtype=np.float))
    meta_data = {"fiducial_wavelength": np.mean(xreflams)}
    print("Avg lam, Std lam: %f, %f" % (np.mean(xreflams),
                                        np.std(xreflams)))

    # make arrays
    xs = np.array(xs, dtype=np.float)
    ys = np.array(ys, dtype=np.float)
    # replace None's with -999
    xs[xs != xs] = -999
    ys[ys != ys] = -999

    # Make a KD-Tree
    tree_dat = np.array([xs, ys], dtype=np.float).T
    tree = KDTree(tree_dat)

    # Get the index of the spaxel closest to the center of the CCD
    ignore, center = tree.query([1024, 1024], 1)

    # Set the coords for the center to 0,0 (origin)
    exts[center].Q_ix = 0
    exts[center].R_ix = 0

    print("IN: %d, EXT: %d, LAM: %d, MDN: %d" % (len(exts), n_tot,
                                                 n_lam, n_mdn))

    def populate_hex(to_populate):
        """ Breadth-first search 

            For each spaxel in the datacube this piece of code identifies
            the relative offset based on the rotation matrix defined
            earlier in the file
        """
        global ncol
        # NOTE: self refers to the extraction with index in
        # the called parameter (to_populate).

        # Query 7 for the object + six surrounding members
        # Generate pixel distances relative to self and indices
        v = np.array([xs[to_populate], ys[to_populate]])
        dists, kix = tree.query(v, 7)
        # Transformation matrix
        tfm = P2H * ROT / np.nanmedian(dists) * 2

        # Trim self reference
        if dists[0] < 2:
            dists = dists[1:]
            kix = kix[1:]

        # Trim to within 70 pixels
        ok = dists < 70
        # dists = dists[ok]
        kix = kix[ok]

        # Get the Q,R of the extraction we are populating
        q_this = exts[to_populate].Q_ix
        r_this = exts[to_populate].R_ix

        # Loop over the current nearest extractions
        for nix in kix:
            # Search around current hex via a recrusive call
            # to populate_hex

            # Current extraction
            nv = np.array([xs[nix], ys[nix]])
            # Offset between self and current
            d = nv-v
            # Get offset in Q, R frame
            dp = np.dot(tfm, d)
            # Integer hex positions so we round
            rnd = np.round(dp)

            # Coincident with center spaxel, so skip
            if rnd[0] == 0 and rnd[1] == 0:
                continue

            # Offset larger than a single hex, so skip
            if np.abs(rnd[0]) > 1 or np.abs(rnd[1]) > 1:
                continue

            # A new hex position
            if exts[nix].Q_ix is None:
                exts[nix].Q_ix = q_this + rnd[0]
                exts[nix].R_ix = r_this + rnd[1]
                # Drill deeper
                populate_hex(nix)
            # We've been to this position before
            else:
                # Check if our hex positions agree
                if (exts[nix].Q_ix != q_this + rnd[0]) or \
                        (exts[nix].R_ix != r_this + rnd[1]):
                    print("collision: %5.1f %5.1f %5.1f %5.1f" %
                          (exts[nix].Q_ix, q_this + rnd[0],
                           exts[nix].R_ix, r_this + rnd[1]))
                    # Update with this one, but don't drill more
                    # otherwise we get in an infinite loop
                    exts[nix].Q_ix = q_this + rnd[0]
                    exts[nix].R_ix = r_this + rnd[1]
                    ncol += 1
        # end loop over current nearest extractions
    # end def populate_hex

    populate_hex(center)

    print("Number of collisions: %d" % ncol)

    # Now convert Q/R to even-Q X/Y

    # Get arrays of hex positions
    qs = np.array([ex.Q_ix for ex in exts], dtype=np.float)
    rs = np.array([ex.R_ix for ex in exts], dtype=np.float)

    # convert to pixel X,Y
    # xs = np.sqrt(3) * (qs + rs/2.0)
    xs = 2.1 * (qs + rs/2.2)
    ys = 1.5 * rs

    # t= np.radians(165.0+45)
    # Hex angle relative to positive Y pixel axis
    t = np.radians(180+13.5)
    # Rotation matrix
    rot = np.array([[np.cos(t), -np.sin(t)],
                    [np.sin(t),  np.cos(t)]])

    # Loop over extractions and project onto X,Y
    for ix, ex in enumerate(exts):
        # Project into X,Y frame
        p = np.dot(rot, np.array([xs[ix], ys[ix]]))
        # Convert to RA, Dec
        # Note 0.633 is plate scale measured on 22 May 2014.
        # hex_scale = 0.315
        ex.Xhex_as = p[0] * 0.300
        ex.Yhex_as = p[1] * 0.355
        ex.X_as = ex.Xhex_as
        ex.Y_as = ex.Yhex_as

    # Write out the cube
    np.save(outname, [exts, meta_data])
    print("Wrote %s.npy" % outname)


def reject_outliers(data, m=2.):
    """Reject outliers beyond `m` sigma

    Args:
        data (numpy float array): data values
        m (float): sigma factor for rejections threshhold

    Returns:
        (numpy float array): input array with outliers removed

    """
    d = np.abs(data - np.nanmedian(data))
    mdev = np.nanmedian(d)
    d[d != d] = 10000.
    s = d/mdev if mdev else 0.
    return data[s < m]


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description="""Convert an extracted file into a data cube.
STEP is either:
make: To create the data cube (once per night)
extract: To extract the cube (one for each observation)
        """, formatter_class=argparse.RawTextHelpFormatter)

    parser.add_argument('extracted', type=str, help='Extracted file')
    parser.add_argument('--step', type=str, default='make',
                        help='[make|extract|dump]')
    parser.add_argument('--outname', type=str, default='cube',
                        help='Output cube name')

    args = parser.parse_args()

    if args.outname is not None:
        args.outname = os.path.splitext(args.outname)[0]

    step = args.step
    infile = args.extracted

    if step == 'make':
        print("\nMAKING cube from %s " % infile)
        ext = np.load(infile)
        extraction_to_cube(ext, outname=args.outname)
    elif step == 'extract':
        print("\nEXTRACTING from %s " % infile)
        ext, meta = np.load(infile)
        qr_to_img(ext, size=2, outname=args.outname)
    elif step == 'dump':
        print("\nDUMPING from %s to %s_dump.txt" % (infile, infile))
        cube = np.load(infile)
        Xs = np.array([c.X_as for c in cube])
        Ys = np.array([c.Y_as for c in cube])
        Sid = np.array([c.seg_id for c in cube])

        dat = np.array([Xs, Ys, Sid])
        np.savetxt("%s_dump.txt" % infile, dat.T)
    else:
        print("NO STEP TO PERFORM")
