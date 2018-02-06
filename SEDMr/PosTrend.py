"""Plot efficiency trend for SEDM

"""
import glob
import os

import pylab as pl
import numpy as np
from astropy.time import Time

sdir = '/scr2/sedmdrp/redux'
fspec = os.path.join(sdir, '20??????')
dlist = sorted([d for d in glob.glob(fspec) if os.path.isdir(d)])[1:]

jd0 = 2457000.0
jd = []
pjd = []
xas = []
yas = []
ras = []
for d in dlist:
    print(d)
    ddate = d.split('/')[-1]
    dtime = Time(ddate[0:4]+'-'+ddate[4:6]+'-'+ddate[6:])

    fspec = os.path.join(d, 'sp_STD-*_obs1*.npy')
    slist = glob.glob(fspec)

    for s in slist:
        ss = np.load(s)[0]
        if 'position' not in ss:
            continue
        if ss['quality'] != 0:
            continue

        xas.append(ss['position'][0])
        yas.append(ss['position'][1])
        ras.append(np.sqrt(ss['position'][0]**2 + ss['position'][1]**2))

        jd.append(dtime.jd)
        pjd.append(dtime.jd - jd0)

pl.plot(pjd, xas, 'x', markersize=2.0, linestyle='None', label='X (RA)')
pl.plot(pjd, yas, 'D', markersize=2.0, linestyle='None', label='Y (Dec)')
pl.plot(pjd, ras, '^', markersize=2.0, linestyle='None', label='R')
pl.xlabel('JD - 2457000')
pl.ylabel('Offset (asec)')
pl.title('Position Trend')
pl.legend(loc=2)
pl.grid(True)
ofil = os.path.join(sdir, 'SEDM_pos_trend.pdf')
pl.savefig(ofil)


