import os
import astropy.io.fits as pf
import NPK.Util as UU
import numpy as np
import warnings

from scipy.interpolate import interp1d
from numpy.polynomial.chebyshev import chebval

# These two are default values needed by WCS to achieve a constant R~100
# wavelength grid for SEDM. 
CRVAL1 = 239.5
CRPIX1 = 88.98


def readspec(path, corrname='std-correction.npy'):
    """Read numpy spec file

    Args:
        path (str): full filespec for numpy sp_*.npy file
        corrname (str): flux correction file name
    
    Returns:
        wavelength array [N]: in nm
        spectrum array [N]: in erg/s/cm2/ang
        sky spectrum array [N]: in spectrum units
        standard deviation of spectrum [N]: in spectrum units
        Spectrum object: full spectrum from where above derived
        meta {}: The meta dictionary associated with the spectrum
        
    """

    # Check for local version
    if not os.path.isfile(corrname):
        # Check SEDM_REF env var
        sref = os.getenv("SEDM_REF")
        if sref is not None:
            corrname = os.path.join(sref, 'std-correction.npy')
        else:
            corrname = '/scr2/sedm/ref/std-correction.npy'

    print "Attempting to load standard correction in: %s" % corrname
        
    ss = np.load(path)[0]

    corr = np.load(corrname)[0]
    corf = interp1d(corr['nm'], corr['correction'], bounds_error=False,
                    fill_value=1.0)

    if 'maxnm' in corr:
        maxnm = corr['maxnm']
    else:
        maxnm = 920.0
    
    lam, spec = ss['nm'], ss['ph_10m_nm']*corf(ss['nm'])

    if 'skyph' in ss:
        skyspec = ss['skyph'] * corf(lam)
    else:
        skyspec = None
        print "Spectrum in %s has no sky spectrum" % path 
    
    if 'var' in ss:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", category=RuntimeWarning)
            std = np.sqrt(np.abs(ss['var']) * corf(lam)*corf(lam))
    else:
        std = None
        print "Spectrum in %s has no variance spectrum" % path

    if 'meta' in ss:
        meta = ss['meta']
    else:
        meta = {}

    meta['maxnm'] = maxnm

    return lam, spec, skyspec, std, ss, meta


def readfits(path):
    """Read fits file at path or path.gz"""

    if not os.path.exists(path):
        if os.path.exists("%s.gz" % path):
            path += ".gz"
        else:
            raise Exception("The file at path %s or %s.gz does not exist" %
                            (path, path))

    hdulist = pf.open(path)
    
    return hdulist


def writefits(towrite, fname, no_lossy_compress=False, clobber=False):

    if type(towrite) == pf.PrimaryHDU:
        hlist = pf.HDUList(towrite)
    elif type(towrite) == pf.HDUList:
        hlist = towrite
    else:
        hlist = pf.HDUList(pf.PrimaryHDU(towrite))

    if no_lossy_compress: 
        if '.gz' in fname:
            hlist.writeto(os.path.splitext(fname)[0], overwrite=clobber)
        else:
            hlist.writeto(fname, overwrite=clobber)
        return
    
    if '.gz' in fname:
        n = os.path.splitext(fname)[0]
    else:
        n = fname
    hlist[0].data = UU.floatcompress(hlist[0].data)
    hlist.writeto(fname, overwrite=clobber)
    
    os.system("gzip  %s" % n)
    

def convert_spectra_to_recarray(spectra):
    """Returns an Numpy recarray version of spectra"""

    keys = spectra[0].__dict__.keys()

    to_remove = ["mdn_coeff", "hg_lines", "spec", "specw"]
    map(keys.remove, to_remove)

    types = []
    for key in keys:
        l = np.size(getattr(spectra[0], key))
        if key == 'lamcoeff':
            l = 6
        if key == 'spec':
            l = 265
        if key == 'specw':
            l = 265

        types.append((key, np.float, l))

    to_handle = []
    for spectrum in spectra:
        res = []
        for kix, key in enumerate(keys):
            toadd = getattr(spectrum, key)

            len_key = types[kix][2]

            if toadd is None:
                val = np.zeros(len_key)
            elif len_key != np.size(toadd):
                val = np.zeros(len_key)
                val[0:np.size(toadd)] = toadd
            else:
                val = toadd

            res.append(val)

        to_handle.append(res)

    ra = np.rec.array(to_handle, dtype=types)
    return ra


def exp_fid_wave(CRVAL1=239.5, CRPIX1=88.98):
    """Return a fiducial wavelength grid appropraite for FITS representation
    
    Computation performed with Mathematica
    """
    
    return CRVAL1 * np.exp((np.arange(265)+CRPIX1)/CRVAL1)


def convert_spectra_to_img(spectra, CRVAL1, CRPIX1):
    
    lfid = exp_fid_wave(CRVAL1=CRVAL1, CRPIX1=CRPIX1)
    img = np.zeros((len(spectra), len(lfid)))
    img2 = np.zeros((len(spectra), len(lfid)))
    img[:] = np.nan
    img2[:] = np.nan

    for ix, spectrum in enumerate(spectra):
        spec = spectrum.specw
        if spec is None:
            continue

        img[ix, 0:len(spec)] = spec

        try:
            lam = chebval(np.arange(*spectrum.xrange), spectrum.lamcoeff)
        except:
            continue
        IF = interp1d(lam, spec, bounds_error=False, fill_value=np.nan)

        img2[ix, :] = IF(lfid)

    return img, img2


def write_cube(spectra):
    """Create a FITS file with all spectra written."""

    recarr = convert_spectra_to_recarray(spectra)

    img, img2 = convert_spectra_to_img(spectra, CRVAL1, CRPIX1)

    f1 = pf.PrimaryHDU(img)
    f2 = pf.ImageHDU(np.zeros((10, 10)))
    t3 = pf.BinTableHDU(recarr)
    f4 = pf.ImageHDU(img2)
    f5 = pf.ImageHDU(np.zeros((10, 10)))

    f4.header['CRVAL1'] = CRVAL1
    f4.header['CRPIX1'] = -CRPIX1
    f4.header['CTYPE1'] = 'WAVE-LOG'
    f4.header['CUNIT1'] = 'NM'

    towrite = pf.HDUList([f1, f2, t3, f4, f5])
    towrite.writeto('test.fits', overwrite=True)

