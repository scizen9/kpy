# -*- coding: utf-8 -*-
"""
Created on Wed Aug  3 15:09:52 2016

@author: nblago
"""


import numpy as np
from scipy import interpolate
from matplotlib import pylab as plt

def extract_percentile(mydir = "20160826", myname = "PTF16fgz.npy", plot=True, write=True, maxperc=30):
    '''
    Test to extract only the spaxels with higher signal to  noise.
    '''
    
    def get_header(myfile):
        header = ""
        with open(myfile) as f:
            ls = f.readlines()
            for li in ls:
                if (li.strip().startswith("#")):
                    header += li.replace("#", "")
                
        return header
    
    
    header = get_header("/scr2/sedmdrp/redux/%s/%s"%(mydir, myname.replace(".npy", "_SEDM.txt")))
    
    #Extraction parameters (chosen spaxels)
    E, meta = np.load("/scr2/sedmdrp/redux/%s/%s"%(mydir, myname))
    #Extraction spectra
    Es = np.load("/scr2/sedmdrp/redux/%s/sp_%s"%(mydir, myname))
    #EsA = np.load("/scr2/sedmdrp/redux/%s/sp_A_%s"%(mydir, myname))
    #EsB = np.load("/scr2/sedmdrp/redux/%s/sp_B_%s"%(mydir, myname))
    std = np.load("/scr2/sedmdrp/redux/%s/std-correction.npy"%mydir)
    dome = np.load("/scr2/sedmdrp/redux/%s/dome.npy"%mydir)
    
    
    mask = ~np.isnan(std[0]["correction"])
    waves = std[0]["nm"][mask]
    wmask = (waves>500)*(waves<900)
    
    for percent in np.arange(00,maxperc,10):
     
        Aflux = []
        stdA = []
        for spax in Es[0]['object_spaxel_ids_A']:
            try: 
                l1, f1 = E[spax].get_flambda()
                lmask = (l1>500)*(l1<900)
                Aflux.append(np.nansum(f1[lmask]))
                stdA.append(np.std(f1[mask]))
            except:
                pass
        maskA = np.array(Aflux)>np.percentile(Aflux, percent)
        
        Bflux = []
        stdB = []
        for spax in Es[0]['object_spaxel_ids_B']:
            try: 
                l1, f1 = E[spax].get_flambda()
                lmask = (l1>500)*(l1<900)
                Bflux.append(-1*np.nansum(f1[lmask]))
                stdB.append(np.std(f1[mask]))
        
            except:
                pass
        maskB = np.array(Bflux)>np.percentile(Bflux, percent)
        
        
        
        medstdA = np.median(stdA)
        medstdB = np.median(stdB)
        
        
        Aspec = []
        Bspec = []
        
        for spax in Es[0]['object_spaxel_ids_A'][maskA]:
            try: 
                l1, f1 = E[spax].get_flambda()
                d1, df1 = dome[0][spax].get_flambda()
                mask_outliers = np.repeat(True, len(l1)) #np.abs(f1 - np.median(f1))<8*medstdA
                specInt = interpolate.interp1d(l1[mask_outliers], f1[mask_outliers], kind="linear", bounds_error=False)
        
                #f1 = f1/df1
                plt.plot(waves,specInt(waves), "r-", alpha=0.05)
                Aspec.append(specInt(waves))
            except IOError:
                pass
            
        for spax in Es[0]['object_spaxel_ids_B'][maskB]:  
            try: 
                l2, f2 = E[spax].get_flambda()
                mask_outliers = np.abs(f1 - np.median(f1))<4*medstdB
        
                specInt = interpolate.interp1d(l2[mask_outliers], f2[mask_outliers], kind="linear", bounds_error=False)
        
                d2, df2 = dome[0][spax].get_flambda()
                #f1 = f1/df1
                plt.plot(waves, specInt(waves), "b-", alpha=0.05)
                Bspec.append(specInt(waves))
            except:
                pass
        
          
        Aspec = np.array(Aspec)
        Bspec = np.array(Bspec)
        
        
        
        sp = interpolate.interp1d(waves, std[0]["correction"][mask], kind="linear", bounds_error=False)
        #spec = spec / np.max(spec) * np.max(np.nanmedian(Aspec, axis=0))
        
        correction = sp(waves)
        
        
        specA = interpolate.interp1d(waves, np.nansum(Aspec, axis=0), kind="linear", bounds_error=False)
        specB = interpolate.interp1d(waves, np.nansum(Bspec, axis=0), kind="linear", bounds_error=False)
        specSum =  (specA(waves) - specB(waves) )
        specSum = specSum / np.max(specSum[wmask]) * np.max(np.nanmedian(Aspec, axis=0))
        
        
        spec =  (specA(waves) - specB(waves) ) * (correction) #/np.median(correction))
        spec = spec / np.max(spec[wmask]) * np.max(np.nanmedian(Aspec, axis=0))
        
        writemask = (waves>400)*(waves<950)
        if (write):
        
            plt.savetxt( "/home/nblago/classifications/%s/%d_%s"%(mydir,percent, myname.replace(".npy", ".txt")), np.array([np.array(waves[writemask][::-1])*10, spec[writemask][::-1] ]).T , fmt="%.1f %.4e", header=header)
        
        if(plot):
            plt.plot(waves, np.nanmedian(Aspec, axis=0), "r-", lw=2, label="Median A")
            plt.plot(waves, np.nansum(Aspec, axis=0)*np.max(np.nanmedian(Aspec, axis=0))/np.max(np.nansum(Aspec, axis=0)), "m--", lw=2, ls="--", label="Sum A")
            plt.plot(waves, -1*np.nanmedian(Bspec, axis=0), "g-", lw=2, label="Median B")
            plt.plot(waves, -1*np.nansum(Bspec, axis=0)*np.max(np.nanmedian(Bspec, axis=0))/np.max(np.nansum(Bspec, axis=0)), "b--", lw=2, ls="--", label="Sum B")
            
            
            plt.plot(waves, np.nansum(Aspec, axis=0)*np.max(np.nanmedian(Aspec, axis=0))/np.max(np.nansum(Aspec, axis=0)) - np.nanmedian(Aspec, axis=0), lw=1, label="Sum - Median A")
            plt.plot(waves, -1*np.nansum(Bspec, axis=0)*np.max(np.nanmedian(Bspec, axis=0))/np.max(np.nansum(Bspec, axis=0)) + np.nanmedian(Bspec, axis=0), lw=1, label="Sum - Median B")
            
        
            plt.plot(waves[writemask], spec[writemask], color="orange", lw=2)
            plt.plot(waves[writemask], specSum[writemask], color="cyan", lw=2)
            
            plt.legend(loc="best", frameon=False)
        
            plt.show()


