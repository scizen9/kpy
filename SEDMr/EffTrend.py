"""Plot efficiency trend for SEDM

"""
import glob
import os

import pylab as pl
import numpy as np
from astropy.time import Time

sdir = '/scr2/sedmdrp/redux'
fspec = os.path.join(sdir, '20??????')
dlist = sorted([d for d in glob.glob(fspec) if os.path.isdir(d)])

jd = []
ef1 = []
ef2 = []
ef3 = []
ef4 = []
ef5 = []
for d in dlist:
    print d
    ddate = d.split('/')[-1]
    dtime = Time(ddate[0:4]+'-'+ddate[4:6]+'-'+ddate[6:])

    fspec = os.path.join(d, 'sp_STD*.npy')
    slist = glob.glob(fspec)

    for s in slist:
        ss = np.load(s)[0]
        if 'efficiency' not in ss:
            continue

        ef = ss['efficiency']*100.
        wl = ss['nm']

        ef1.append(np.nanmean(ef[(wl > 400) * (wl < 500)]))
        ef2.append(np.nanmean(ef[(wl > 500) * (wl < 600)]))
        ef3.append(np.nanmean(ef[(wl > 600) * (wl < 700)]))
        ef4.append(np.nanmean(ef[(wl > 700) * (wl < 800)]))
        ef5.append(np.nanmean(ef[(wl > 800) * (wl < 900)]))

        jd.append(dtime.jd)

pl.plot(jd, ef1, label='400-500 nm')
pl.plot(jd, ef2, label='500-600 nm')
pl.plot(jd, ef3, label='600-700 nm')
pl.plot(jd, ef4, label='700-800 nm')
pl.plot(jd, ef5, label='800-900 nm')
pl.xlabel('JD')
pl.ylabel('Efficiency(%)')
pl.title('Efficiency Trend')
pl.legend(loc=2)
pl.grid(True)
ofil = os.path.join(sdir, 'SEDM_eff_trend.pdf')
pl.savefig(ofil)


