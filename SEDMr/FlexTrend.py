"""Plot flexure trend for SEDM

"""
import glob
import os
import csv
from itertools import izip

import pylab as pl
import numpy as np
from astropy.time import Time

sdir = '/scr2/sedmdrp/redux'
fspec = os.path.join(sdir, '20??????')
dlist = sorted([d for d in glob.glob(fspec) if os.path.isdir(d)])[1:]

jd0 = 2457000.0
jd = []
flxy = []
flxx = []
ffs = []

for d in dlist:
    print(d)
    ddate = d.split('/')[-1]
    dtime = Time(ddate[0:4]+'-'+ddate[4:6]+'-'+ddate[6:])

    fspec = os.path.join(d, 'flex_*.npy')
    slist = glob.glob(fspec)

    for s in slist:
        ss = np.load(s)[0]
        if 'dXnm' not in ss or 'dYpix' not in ss:
            continue

        dx = ss['dXnm']
        dy = ss['dYpix']

        # Fill flexure vectors
        flxx.append(dx)
        flxy.append(dy)

        ffs.append(s)

        jd.append(dtime.jd)

pl.plot(jd-jd0, flxx, 'bD', ms=2.0, label='dX nm')
pl.plot(jd-jd0, flxy, 'g^', ms=2.0, label='dY pix')

pl.xlabel('JD - 2457000')
pl.ylabel('Flexure')
pl.title( 'Flexure Trend')
pl.legend(loc=2)
pl.grid(True)
pl.ylim(-2, 2)
ofil = os.path.join(sdir, 'SEDM_flex_trend.pdf')
pl.savefig(ofil)
with open('SEDM_flex_trend.txt', 'wb') as dfil:
    dat_writer = csv.writer(dfil, delimiter=" ", quoting=csv.QUOTE_MINIMAL)
    dat_writer.writerows(izip(jd, flxx, flxy, ffs))
