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
pjd = []
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
        if 'xfwhm' not in ss or 'yfwhm' not in ss:
            continue
        if ss['xfwhm'] == 0. or ss['yfwhm'] == 0.:
            continue

        dx = ss['dXnm']
        dy = ss['dYpix']

        # Fill flexure vectors
        flxx.append(dx)
        flxy.append(dy)

        ffs.append(s)

        jd.append(dtime.jd)
        pjd.append(dtime.jd-jd0)

pl.plot(pjd, flxx, 'bD', markersize=2.0, linestyle='None', label='dX nm')
pl.plot(pjd, flxy, 'g^', markersize=2.0, linestyle='None', label='dY pix')

pl.xlabel('JD - 2457000')
pl.ylabel('Flexure')
pl.title( 'Flexure Trend')
pl.legend(loc=2)
pl.grid(True)
pl.ylim(-5, 5)
ofil = os.path.join(sdir, 'SEDM_flex_trend.pdf')
tfil = os.path.join(sdir, 'SEDM_flex_trend.txt')
pl.savefig(ofil)
with open(tfil, 'wb') as dfil:
    dat_writer = csv.writer(dfil, delimiter=" ", quoting=csv.QUOTE_MINIMAL)
    dat_writer.writerows(izip(jd, flxx, flxy, ffs))
