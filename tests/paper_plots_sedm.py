# -*- coding: utf-8 -*-
"""
Created on Wed Mar  1 18:38:17 2017

@author: nblago
"""
import numpy as np
from matplotlib import pylab as plt
from numpy.lib import recfunctions as rfn
import os
import sextractor
import spec_utils
import glob

def add_fwhm(myfile, outfile):

    a = np.genfromtxt(myfile, dtype=None, names=True, delimiter=",")
    fwhm = np.zeros(len(a))
    

    for i, f in enumerate(a["filename"]):
        sexfile = f.replace("reduced/", "reduced/sextractor/").replace(".fits", ".sex")

        if (not os.path.isfile(sexfile)):
            sextractor.run_sex([f], mask=False, cosmics=False, overwrite=False)
               
        if (os.path.isfile(sexfile)):
            s = np.genfromtxt(sexfile, comments="#")
            #Select round sources (ellipticity is 1-axis_ratio)
            s = s[s[:,8]<np.percentile(s[:,8], 30)]
            #Select bright magnitudes
            s = s[s[:,4]<np.percentile(s[:,4], 20)]
            
            fwhm[i]= np.percentile(s[:,7], 33)*0.394

            
    out = rfn.append_fields(a, names="fwhm", data=fwhm, usemask=False)
    
    if (not outfile is None):
        np.savetxt(outfile, out, header="#object,filename,filter,std,stderr,inst,insterr,jd,airmass,color,exptime,fwhm", fmt="%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.3f,%.3f,%.2f,%.3f")
        
    return out
    
def compute_S2N(allstarsfile="/home/nblago/Documents/allstars_2017_fwhm.log", exptime=180):
    '''
    Receives as input the file of all stars taken during the observing night.
    '''
    
    a = np.genfromtxt(allstarsfile, dtype=None, names=True, delimiter=",")
    a = a[a["exptime"]==exptime]    
    a = a[~np.isnan(a["exptime"])]
    a = a[~np.isnan(a["std"])]
    a = a[~np.isnan(a["insterr"])]
    a = a[ (a["std"] - a["inst"]) < 24]
    a = a[a["airmass"] < 2.]
    a = a[a["fwhm"]!=0]
    a = a[a["fwhm"]<3]
    
    if "fwhm" not in a.dtype.names:
        a = add_fwhm(a)
    
    if (len(a)==0):
        print "No registers to show"
        return
        
    S = 10**(-0.4 * a["inst"]) * a["exptime"]
    
    N = S - (S / (10**(0.4*a["insterr"])))
    
    
    print zip( a["filename"][(a["std"]<18) * (S/N)<7.5], a["jd"][(a["std"]<18) * (S/N)<7.5])
    print len(a["filename"][(a["std"]<18) * (S/N)<7.5])
    mask_bad = (a["std"]<18) * ((S/N)<7.5)
    
    f =plt.figure()
    
    colors = ['purple', 'green', 'red', 'orange']
    for i, filt in enumerate(['u', 'g', 'r', 'i']):
        mask_filter = a['filter']==filt
        ax1 = plt.subplot2grid((2,2), (i/2, i%2), colspan=1)
        scat = ax1.scatter(a["std"][mask_filter*~mask_bad], (S/N)[mask_filter*~mask_bad], s=2, c=a["fwhm"][mask_filter*~mask_bad])
        ax1.text(0.85, 0.8, filt, transform=ax1.transAxes, fontsize=16)
        ax1.set_ylim(1,700)
        ax1.set_xlim(15,21.1)
        ax1.set_yscale("log")
        if (i%2==1):
            ax1.set_yticklabels("")
        if (i%2==0):
            ax1.set_ylabel("log$_{10}$ (S/N)")
        if ( i/2 == 1):
            ax1.set_xlabel("mag")
        if ( i/2 == 0):
            ax1.set_xticklabels("")
        plt.legend()    
    
    
    plt.subplots_adjust(wspace=0.0, hspace=0, right=0.85)
    cbar_ax = f.add_axes([0.87, 0.15, 0.02, 0.7])
    cb = f.colorbar(scat, cax=cbar_ax, boundaries=[1.0, 1.2, 1.4, 1.6, 1.8, 2 ,2.2 , 2.4, 2.6, 2.8, 3])
    cb.set_label("fwhm")
    
    plt.savefig("/home/nblago/Projects/SEDM/plots_paper/signal2noise_%s.pdf"%myexptime, dpi=200)
    plt.show()
    
    '''mask_sn = ~np.isnan(S/N)
    plt.figure()
    plt.hist2d(a["std"][~mask_bad*mask_sn], (S/N)[~mask_bad*mask_sn], range=[[16,22], [0, 200]], bins=30)
    plt.show()'''
    
def compute_S2N_mags(allstarsfile="/home/nblago/Documents/allstars_2017_fwhm.log", exptime=180):
    '''
    Receives as input the file of all stars taken during the observing night.
    '''
    
    a = np.genfromtxt(allstarsfile, dtype=None, names=True, delimiter=",")
    a = a[a["exptime"]==exptime]    
    a = a[~np.isnan(a["exptime"])]
    a = a[~np.isnan(a["std"])]
    a = a[~np.isnan(a["insterr"])]
    a = a[ (a["std"] - a["inst"]) < 24]
    a = a[a["airmass"] < 2.]
    a = a[a["fwhm"]!=0]
    a = a[(a["fwhm"]<2.5)*(a["fwhm"]>1.5)]
    
    if "fwhm" not in a.dtype.names:
        a = add_fwhm(a)
    
    if (len(a)==0):
        print "No registers to show"
        return
        
    S = 10**(-0.4 * a["inst"]) * a["exptime"]
    
    N = S - (S / (10**(0.4*a["insterr"])))
    
    
    print zip( a["filename"][(a["std"]<18) * (S/N)<7.5], a["jd"][(a["std"]<18) * (S/N)<7.5])
    print len(a["filename"][(a["std"]<18) * (S/N)<7.5])
    mask_bad = (a["std"]<18) * ((S/N)<7.5)
    
    f =plt.figure(figsize=(12,8))
    
    colors = ['green', 'red']
    widths = [4,2.5,1]
    for i, filt in enumerate(['g', 'r']):
        for j, m in enumerate([17.5, 18.5, 19.5]):
            ax = plt.subplot2grid((2, 3), (i, j), colspan=1)
            mask_filter_mag = (a['filter']==filt) * (np.abs(a['std']-m)<0.5)
            plt.hist((S/N)[mask_filter_mag*~mask_bad], bins=15, color=colors[i], histtype="step", label="%.1f"%m, range=(0,200), normed=True, linewidth=widths[j]) 
            ax.set_xlabel("S/N")
            if (j==0):
                ax.set_ylabel("Frequency")
            else:
                ax.set_yticklabels("")
            plt.legend()
    
    
    plt.subplots_adjust(wspace=0.05)   
    #plt.subplots_adjust(wspace=0.1, hspace=0.1, right=0.85)   
    plt.savefig("/home/nblago/Projects/SEDM/plots_paper/signal2noise_%s_gr.pdf"%(exptime), dpi=200)
    plt.show()
    
    '''mask_sn = ~np.isnan(S/N)
    plt.figure()
    plt.hist2d(a["std"][~mask_bad*mask_sn], (S/N)[~mask_bad*mask_sn], range=[[16,22], [0, 200]], bins=30)
    plt.show()'''
    
def plot_two_spec(folder="/home/nblago/Projects/SEDM/plots_paper/specs/"):

    
    plt.figure(figsize=(6,8))

    files = glob.glob(os.path.join(folder, "*ascii"))
    files.sort()
    
    for i, f in enumerate(files):
        s1 = spec_utils.load_spec(f)
        exptime = 0
        with open(f, "r") as rf:
            rlines = rf.readlines()[0:20]
            for l in rlines:
                if ("EXPTIME" in l):
                    exptime = float(l.split(":")[1])
            
        base = os.path.basename(f)
        n, name, date, clas, z = (base.replace(".ascii", "")).split("_")
        z = float(z)
        s1[:,0] = s1[:,0]/(1.+z)
        mask = (s1[:,0]<9000)*(s1[:,0]>4000)
        s1 = s1[mask]
        
        y = s1[:,1] #spec_utils.smooth(s1[:,0], s1[:,1], 1)
        y = y - np.median(y)
        y = y/(np.max(y) - np.min(y))
        
        plt.step(s1[:,0], y + i*0.75+0.75, lw=2, label="%s %ds"%(clas, exptime))
        
        #leg = plt.legend(loc="lower left", framealpha=1., frameon=True)
        #leg.get_frame().set_linewidth(0.0)
    
    ymin, ymax = plt.ylim()
    plt.ylim(ymin=0, ymax=ymax*1.1)
    plt.xlim(xmin=3900)
    plt.legend(loc="lower right")
    spec_utils.plot_element_lines("/home/nblago/workspace/datared/data/spectra/ptf16fnl_lines_simple.txt", ymin, ymax*1.1)   

    plt.gca().minorticks_on()
    plt.xlabel("Restframe wavelength [$\AA$]")
    plt.ylabel("Flux F$_{\lambda}$(erg/$\AA$/cm$^2$/s) + constant")
    plt.tight_layout()
    plt.savefig("/home/nblago/Projects/SEDM/plots_paper/spectra_sedm.pdf", ddp=200)
    plt.show()