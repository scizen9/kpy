import numpy as np
from astropy.coordinates import Angle

"""
Struct:
http://stackoverflow.com/questions/1305532/convert-python-dict-to-object
"""


class Struct(object):
    """Converts a dictionary into an object

    Example:
        a = {'a': 1, 'b': 2}
        o = Struct(a)
        print o.a, o.b

    """
    def __init__(self, **entries):
        """s = Struct(***dict)"""
        self.__dict__.update(entries)


def floatcompress(data, ndig=14):
    """Adapted from Finkbeiner IDL routine floatcompress"""

    t = data.dtype
    if not ((t == 'float32') or (t == 'float64')):
        raise Exception("Only works on floating point numbers")

    wzer = np.where(data == 0)
    data[wzer] = 1.0

    log2 = np.ceil(np.log(np.abs(data)) / np.log(2.0))
    mant = np.round(data/2.0**(log2 - ndig))/2.0**ndig
    out = mant*2.0**log2

    out[wzer] = 0.0
    return out


def par_angle(ha, dec, lat=Angle("33d 21m 21.6s")):
    """ Compute parallactic angle given an hour angle and declination.

    Parallactic angle is defined such that 180 is N/S.

    Latitude defaults to Palomar. Must use Astropy Angle class"""

    if type(ha) != Angle:
        ha = Angle(ha, unit='hour')
    if type(dec) != Angle:
        dec = Angle(dec, unit='deg')

    return Angle(np.arctan2(-np.sin(ha),
                            np.cos(dec) * np.tan(lat) - np.sin(dec) *
                            np.cos(ha)))


def air_index(lam, p=600, t=7, f=8):
    """ Returns index of refraction of air-1 at
        lam in micron at vacuum
        p is pressure in mm Hg
        t is temperature in deg C
        f is water vapor pressure in mm Hg
    """

    k1 = (1/lam)**2
    nm1e6 = 64.328 + 29498.1/(146-k1) + 255.4/(41-k1)

    nm1e6 *= p * (1 + (1.049-0.0157 * t)*1e-6*p) / (720.883 * (1 +
                                                               0.003661 * t))

    nm1e6 -= 0.0624 - 0.000680 * k1 / (1 + 0.003661 * t) * f

    return nm1e6/1e6


def atm_disper(l2, l1, airmass, **kwargs):
    """ atmospheric dispersion in arcsecond between l2 and l1 in micron
        at a given airmass. See air index for documentation on pressure,
        temperature, and water vapor pressure"""

    z = np.arccos(1.0/airmass)
    return 206265 * (air_index(l2, **kwargs) - air_index(l1,
                                                         **kwargs)) * np.tan(z)
