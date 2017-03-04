"""Spectra class definition

Functions
    * :func:`find_ha`   return index closest to Halpha

"""
import numpy as np
import scipy.spatial
from numpy.polynomial.chebyshev import chebval
import warnings


def find_ha(cc):
    """Find out where in X position Halpha falls

    Args:
        cc (float vector): coefficients of fit giving
            wavelength versus spatial position

    Returns:
        float: spatial position corresponding to Halpha

    Note:
        assumes spatial position falls between 30 and 100
        arcsec and has an accuracy limit of 0.1 arcsec

    """
    ix = np.arange(30, 100, .1)

    haix = np.argmin(np.abs(chebval(ix, cc) - 656.3))

    return ix[haix]


class Spectra(object):
    """Spectra class"""
    KT = None       # The KD Tree
    data = None     # The actual data
    good_positions = []     # The mapping of KT data to data so that
                            # some ix in the KT.data is the same
                            # as data[KT.good_positions[ix]]

    def __init__(self, data=None):

        positions = []
        good_positions = []
        self.data = data

        for ix, el in enumerate(data):

            try:
                l, fl = el.get_counts()
                fl = el.specw
                haix = find_ha(el.lamcoeff)
            except:
                continue

            good_positions.append(ix)

            x = el.X_as
            y = el.Y_as

            positions.append((x, y))

        if len(positions) == 0:
            raise Exception("For some reason, no good spectrum exists"
                            " in the submitted spectral set.")

        data = np.array(positions)
        bad = (data != data)
        data[bad] = -999
        self.KT = scipy.spatial.KDTree(data)
        self.good_positions = np.array(good_positions)

    def to_xyv(self, lmin=500., lmax=700.):
        """ Method converts a spectral set into X, Y, and intensity tuples
        
        X
            x position relative to center of IFU in arcsec
            based on the segmentation map.

        Y
            y position relative to center of IFU in arcsec
            based on the segmentation map.

        V
            median signal strength in the lmin to lmax range.

        Args:
            lmin (float): minimum wavelength in nanometers
            lmax (float): maximum wavelength in nanometers

        Returns:
            X,Y,V tuple: the X location (in arcsec), Y location
            (in arcsec), and the median value (V) of the spaxel
            between lmin and lmax

        """

        Xs = []
        Ys = []
        Vs = []

        for ix, XY in enumerate(self.KT.data):

            datix = self.good_positions[ix]
            el = self.data[datix]
            try:
                l, fl = el.get_flambda()
            except: 
                continue

            ok = (l > lmin) & (l <= lmax)

            Xs.append(XY[0])
            Ys.append(XY[1])
            with warnings.catch_warnings():
                warnings.simplefilter("ignore", category=RuntimeWarning)
                Vs.append(np.nanmedian(el.spec[ok]))

        return (np.array(Xs),
                np.array(Ys),
                np.array(Vs))
