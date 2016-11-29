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
import pyfits as pf
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


def QR_to_img(exts, Size=4, outname="cube.fits"):
    """Convert a data cube to a fits image

    Args:
        exts (list of Extraction): extractions to convert (see Extraction.py)
        Size (int): expansion factor, defaults to 4
        outname (str): output fits file name, defaults to cube.fits

    Returns:
        None

    """
    
    Xs = np.array([ext.X_as for ext in exts], dtype=np.float)
    Ys = np.array([ext.Y_as for ext in exts], dtype=np.float)

    minx = Size * np.nanmin(Xs)
    miny = Size * np.nanmin(Ys)
    maxx = Size * np.nanmax(Xs)
    maxy = Size * np.nanmax(Ys)

    Dx = (maxx-minx)/.25
    Dy = (maxy-miny)/.25
    l_grid = Wavelength.fiducial_spectrum()
    l_grid = l_grid[::-1]
    dl_grid = np.diff(l_grid)
    l_grid = l_grid[1:]

    img = np.zeros((Dx, Dy, len(l_grid)/2))
    img[:] = np.nan

    XSz = img.shape[0]/2
    YSz = img.shape[1]/2

    allspec = np.zeros((len(exts), len(l_grid)/2))
    for cnt, ext in enumerate(exts):
        if ext.xrange is None:
            continue
        if ext.exptime is None:
            ext.exptime = 1
        if ext.lamcoeff is None:
            continue

        ix = np.arange(*ext.xrange)
        l = chebval(ix, ext.lamcoeff)
        s = ext.specw

        f = interp1d(l, s, fill_value=np.nan, bounds_error=False)
        fi = f(l_grid) / dl_grid
        fi = fi[0:len(fi):2] + fi[1:len(fi):2]

        allspec[cnt, :] = fi

        x = (ext.X_as - minx)/0.25
        y = (ext.Y_as - miny)/0.25

        try: 
            for dx in [-1, 0, 1]:
                for dy in [-1, 0, 1]:
                    img[x+dx, y+dy, :] = fi
        except:
            pass

        outstr = "\rX = %+10.5f, Y = %+10.5f" % (x, y)
        print outstr,
        sys.stdout.flush()

    back = np.median(allspec, 0)

    if 'fits' not in outname:
        outname += '.fits'

    ff = pf.PrimaryHDU(img.T)
    ff.writeto(outname)

    for cnt, ext in enumerate(exts):
        if ext.xrange is None:
            continue
        if ext.exptime is None:
            ext.exptime = 1
        if ext.lamcoeff is None:
            continue

        ix = np.arange(*ext.xrange)
        l = chebval(ix, ext.lamcoeff)
        s = ext.specw 

        f = interp1d(l, s, fill_value=np.nan, bounds_error=False)
        fi = f(l_grid)/dl_grid 
        fi = fi[0:len(fi):2] + fi[1:len(fi):2] - back

        try:
            x = np.round(Size*np.sqrt(3.) * (ext.Q_ix + ext.R_ix/2)) + XSz
            y = np.round(Size*3./2. * ext.R_ix) + YSz
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
        * Z_as: Z position in arcsec (the Z coordinate runs 45 degrees to X and is not a `3rd` dimension).
        * Q_ix: The axial Q coordinate in integral units
        * R_ix: The axial R coordinate in integral units

        The relationship of Q/R to X/Y is defined through the pixel mapping
        matrix times the rotation matrix:

        * P2H = np.array([[np.sqrt(3)/3, -1/3.], [0, 2/3.]])
        * Rot (22 degree)

    """

    global ncol

    # start with our x,y pixel positions
    Xs = [None] * len(exts)
    Ys = [None] * len(exts)

    # these are the hex axial coords
    for ext in exts:
        ext.Q_ix = None
        ext.R_ix = None
    
    # counters for which wavelength solution is used:
    # mdn = median, lam = lambda, tot = total used
    n_mdn = 0
    n_lam = 0
    n_tot = 0

    for ix, ext in enumerate(exts):
        # The X/Y location of a lenslet is based on its
        # trace Y position and where the fiducial wavelength
        # is expected in the X direction.

        # skip bad solutions
        if not ext.ok: continue

        # initialize
        Xs[ix] = -999
        Ys[ix] = -999

        # lamcoeff has precedence
        if ext.lamcoeff is not None:
            coeff = ext.lamcoeff
            n_lam += 1
        elif ext.mdn_coeff is not None: 
            coeff = ext.mdn_coeff
            n_mdn += 1
        else:
            continue

        n_tot += 1

        # get pixel and wavelength vectors
        ixs = np.arange(*ext.xrange)
        LL = chebval(ixs, coeff)
        
        # fill in the x,y pix positions
        Xs[ix] = np.nanmin(ext.xrange) + ext.xrefpix
        Ys[ix] = np.nanmean(ext.yrange)
        ext.X_pix = Xs[ix]
        ext.Y_pix = Ys[ix]
        if ext.xrefpix is not None:
            if ext.xrefpix < len(LL):
                ext.xreflam = LL[ext.xrefpix]
    # End loop over all extractions

    # Record wavelength mean and report stats
    xreflams = reject_outliers(
            np.array([ext.xreflam for ext in exts], dtype=np.float))
    meta = {"fiducial_wavelength": np.mean(xreflams)}
    print "Avg lam, Std lam: %f, %f" % (np.mean(xreflams), np.std(xreflams))

    # make arrays
    Xs = np.array(Xs, dtype=np.float)
    Ys = np.array(Ys, dtype=np.float)
    # replace None's with -999
    Xs[Xs != Xs] = -999
    Ys[Ys != Ys] = -999

    # Make a KD-Tree
    dat = np.array([Xs, Ys], dtype=np.float).T
    tree = KDTree(dat)

    # Get the index of the spaxel closest to the center of the CCD
    ignore, Center = tree.query([1024, 1024], 1)

    # Set the coords for the center to 0,0 (origin)
    exts[Center].Q_ix = 0
    exts[Center].R_ix = 0

    print "IN: %d, EXT: %d, LAM: %d, MDN: %d" % (len(exts), n_tot, n_lam, n_mdn)

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
        v = np.array([Xs[to_populate], Ys[to_populate]])
        Dists, Ixs = tree.query(v, 7)
        # Transformation matrix
        Tfm = P2H * ROT / np.median(Dists) * 2

        # Trim self reference
        if Dists[0] < 2: 
            Dists = Dists[1:]
            Ixs = Ixs[1:]

        # Trim to within 70 pixels
        ok = Dists < 70
        Dists = Dists[ok]
        Ixs = Ixs[ok]

        # Get the Q,R of the extraction we are populating
        q_this = exts[to_populate].Q_ix
        r_this = exts[to_populate].R_ix

        # Loop over the current nearest extractions
        for nix in Ixs:
            # Search around current hex via a recrusive call
            # to populate_hex

            # Current extraction
            nv = np.array([Xs[nix], Ys[nix]])
            # Offset between self and current
            D = nv-v
            # Get offset in Q, R frame
            dp = np.dot(Tfm, D)
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
                    print "collision: ",
                    print exts[nix].Q_ix, q_this + rnd[0], " ",
                    print exts[nix].R_ix, r_this + rnd[1]
                    # Update with this one, but don't drill more
                    # otherwise we get in an infinite loop
                    exts[nix].Q_ix = q_this + rnd[0]
                    exts[nix].R_ix = r_this + rnd[1]
                    ncol += 1
        # end loop over current nearest extractions
    # end def populate_hex

    populate_hex(Center)

    print "Number of collisions: %d" % ncol

    # Now convert Q/R to even-Q X/Y

    # Get arrays of hex positions
    Qs = np.array([ext.Q_ix for ext in exts], dtype=np.float)
    Rs = np.array([ext.R_ix for ext in exts], dtype=np.float)

    # convert to pixel X,Y
    Xs = np.sqrt(3) * (Qs + Rs/2.0)
    Ys = 3/2 * Rs

    # t= np.radians(165.0+45)
    # Hex angle relative to positive Y pixel axis
    t = np.radians(180+22)
    # Rotation matrix
    Rot = np.array([[np.cos(t), -np.sin(t)],
                    [np.sin(t),  np.cos(t)]])

    # Loop over extractions and project onto X,Y
    for ix, ext in enumerate(exts):
        # Project into X,Y frame
        p = np.dot(Rot , np.array([Xs[ix], Ys[ix]]))
        # Convert to RA, Dec
        # Note 0.633 is plate scale measured on 22 May 2014.
        ext.X_as = p[0] * -0.633
        ext.Y_as = p[1] * 0.633

    # Write out the cube
    np.save(outname, [exts, meta])
    print "Wrote %s.npy" % outname


def reject_outliers(data, m=2.):
    """Reject outliers beyond `m` sigma

    Args:
        data (numpy float array): data values
        m (float): sigma factor for rejections threshhold

    Returns:
        (numpy float array): input array with outliers removed

    """
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
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
    parser.add_argument('--outname', type=str, help='Output cube name')

    args = parser.parse_args()

    if args.outname is not None:
        args.outname = os.path.splitext(args.outname)[0]

    step = args.step
    infile = args.extracted

    if step == 'make':
        print "\nMAKING cube from %s " % infile
        ext = np.load(infile)
        cube = extraction_to_cube(ext, outname=args.outname)
    elif step == 'extract':
        print "\nEXTRACTING from %s " % infile
        ext,meta = np.load(infile)
        QR_to_img(ext, Size=2, outname=args.outname)
    elif step == 'dump':
        print "\nDUMPING from %s to %s_dump.txt" % (infile, infile)
        cube = np.load(infile)
        Xs = np.array([c.X_as for c in cube])
        Ys = np.array([c.Y_as for c in cube])
        Sid = np.array([c.seg_id for c in cube])

        dat = np.array([Xs,Ys,Sid])
        np.savetxt("%s_dump.txt" % infile, dat.T)
    else:
        print "NO STEP TO PERFORM"

