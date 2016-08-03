"""Loads standard spectra from Oke 1990 (Aj, 99, 1621).

See:
    ftp://ftp.eso.org/pub/stecf/standards/okestan/aaareadme.oke

Defines:
    Standards[standard_name]: [Wavelength, flux, flux, bin]
    units: 'Wavelength [A], Flux [erg/s/cm/cm/A 10**16], flux [mJy], bin [A]'
"""

import os
import numpy as np


units = 'Wavelength [A], Flux [erg/s/cm/cm/A 10**16], flux [mJy], bin [A]'

sdir = os.path.join(os.path.dirname(__file__), 'standard_stars')

files = os.listdir(sdir)

Standards = {}

for fil in files:
    if fil[0] != 'f':
        continue

    std_name = fil[1:-4]
    dat = np.loadtxt(os.path.join(sdir, fil))
    Standards[std_name] = dat

