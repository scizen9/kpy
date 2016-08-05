# -*- coding: utf-8 -*-
"""
Created on Mon May  9 12:25:51 2016

@author: nblago
"""

import glob
from pyraf import iraf
import numpy as np
import fitsutils
import pyfits as pf
from matplotlib import pylab as plt


tests = {}
for f in glob.glob("ifu*fits"):
    testname = fitsutils.get_par(f, "OBJECT").split()[8]
    if (not tests.has_key(testname)):
        tests[testname] = []
    tests[testname].append(f)
    
names = []
for i, k  in enumerate(tests.keys()):
    name = "comb%d"%i
    np.savetxt(name, tests[k], fmt="%s")
    names.append(name)

for i in glob.glob("comb*"):
    iraf.imcombine("@"+i, "im%s"%i)

images = glob.glob("imcomb[0-9].fits")
images.extend(glob.glob("imcomb[0-9][0-9].fits"))

elevations = [fitsutils.get_par(im, "TEL_EL") for im in images]
azimuths = [fitsutils.get_par(im, "TEL_AZ") for im in images]
jds = [fitsutils.get_par(im, "JD") for im in images]

zipped =zip(images, elevations, azimuths,jds )
zipped.sort(key = lambda t: t[3])

images = [z[0] for z in zipped]

for i in np.arange(len(images)):
    for j in np.arange(i,len(images)):
        # Running IRAF
        iraf.noao(_doprint=0)
        iraf.imred(_doprint=0)
        iraf.ccdred(_doprint=0)
        out = "%s-%s"%(images[i],images[j])
        print images[i], images[j]
        iraf.imarith(images[i], "/", images[j], out)

        

def plot_position(X,Y):

    plt.suptitle("Flexure test %d, %d")
    plt.figure(figsize=(30,30))
    
    dimX = len(images)
    for i in np.arange(len(images)):
        for j in np.arange(i,len(images)):
    
            out = "%s-%s"%(images[i],images[j])
            ax = plt.subplot2grid((dimX,dimX),(i, j))        
            data = pf.open(out)[0].data
            data = data[X-25:X+25, Y-25:Y+25]
            ax.imshow(data, origin="lower", cmap=plt.get_cmap("bwr"), interpolation="none", vmin=0.90, vmax=1.1)#vmin=np.percentile(data, 5), vmax=np.percentile(data,95))
            if(i==0):
                ax.set_title("%.1f, %.1f"%(zipped[j][2], zipped[j][1]) )
            if(i==j):
                ax.set_ylabel("%.1f, %.1f"%(zipped[i][2], zipped[i][1]) )
                ax.set_xlabel("%.1f, %.1f"%(zipped[j][2], zipped[j][1]) )
            
    
    plt.savefig("Flexure_test_%d_%d.png"%(X,Y))
    plt.close("all")
    
plot_position(550,550)
plot_position(1750,550)
plot_position(550,1750)
plot_position(1750,1750)
plot_position(1025,1025)


