# -*- coding: utf-8 -*-
"""
Created on Thu Jun  8 15:50:05 2017

@author: nblago
"""
import numpy as np
import scipy
from scipy.interpolate import interp1d
import astropy.constants as cnt
from astropy import units as u
from matplotlib import pylab as plt
from astropy import stats

class Throughput():
    
    def __init__(self, uimage=None):
        self.uimage = None
        self.gimage = None
        self.rimage = None
        self.iimage = None
        
        self.palomarext = ""
        self.stds = {}
        self.A = 0.8 * np.pi /4 * (1.5e2**2)  #TO BE CHECKED (unobscrued fraction)
    
        self.filters = np.genfromtxt("/home/nblago/Projects/SEDM/scripts/AstrodonSloanGen2Transmission.csv", delimiter=",", names=True)
        self.qe = np.genfromtxt("/home/nblago/Projects/SEDM/scripts/QE.txt", delimiter="\t", names=True)
        self.qesmooth = scipy.interpolate.UnivariateSpline( self.qe["wl"], self.qe["QE55C"] )
        self.qesmooth.set_smoothing_factor(0.5)
        
        # Palomar Extinction Data from Hayes & Latham 1975
        # (Wavelength in Angstroms, Magnitudes per airmass)
        self.palextinct = [
            (3200, 1.058),
            (3250, 0.911),
            (3300, 0.826),
            (3350, 0.757),
            (3390, 0.719),
            (3448, 0.663),
            (3509, 0.617),
            (3571, 0.575),
            (3636, 0.537),
            (3704, 0.500),
            (3862, 0.428),
            (4036, 0.364),
            (4167, 0.325),
            (4255, 0.302),
            (4464, 0.256),
            (4566, 0.238),
            (4785, 0.206),
            (5000, 0.183),
            (5263, 0.164),
            (5556, 0.151),
            (5840, 0.140),
            (6055, 0.133),
            (6435, 0.104),
            (6790, 0.084),
            (7100, 0.071),
            (7550, 0.061),
            (7780, 0.055),
            (8090, 0.051),
            (8370, 0.048),
            (8708, 0.044),
            (9832, 0.036),
            (10255, 0.034),
            (10610, 0.032),
            (10795, 0.032),
            (10870, 0.031)]
        
        self.palextinct = np.array(self.palextinct)
        self.ext = interp1d(self.palextinct[:, 0], self.palextinct[:, 1], kind='cubic', bounds_error=False)

        self.lameff = {"u" : 3560, "g":4830, "r":6260, "i":7670}
        self.deltalam = {"u" : 463, "g": 988, "r":955, "i":1064}

        #self.lameff = {"u" : 3515, "g":4745, "r":6280, "i":8101}
        #self.deltalam = {"u" : 610.0, "g": 1450, "r":1300, "i":1490}
        
        self.lameff = {"u" : 3557, "g":4702, "r":6175, "i":7491}
        self.deltalam = {"u" : 463, "g": 988, "r":955, "i":1064}
        
        '''self.lameff = {}
        self.deltalam = {}
        
        for f in ["u", "g", "r", "i"]:
            mask = self.filters[f]>90
            self.lameff[f] = np.median(self.filters[mask]["Wavelength_nm"])*10
            self.deltalam[f] = (np.max(self.filters[mask]["Wavelength_nm"]) - np.min(self.filters[mask]["Wavelength_nm"])) * 10'''
        
    def get_Palomar_extinction(self, wave, airmass):
        #Read extinction in Palomar. Multiply by airmass
        return self.ext(wave)*airmass
        
        
    def get_avg_Palomar_extinction(self, f, airmass):
        #Read extinction in Palomar. Multiply by airmass
        winit = self.lameff[f] - self.deltalam[f]/2
        wend = self.lameff[f] + self.deltalam[f]/2
        
        wl = np.linspace(winit, wend)
        
        return np.median(self.ext(wl))*airmass #trapz(self.ext(wl)*airmass, wl)
        
    def get_QE_filter(self, f):
        '''
        Returns the average QE in that filter.
        '''
    
        if (f!="i"):
            mask = (self.filters[f] >0.05)*(self.filters["Wavelength_nm"]<864)
        else:
            mask = (self.filters[f] >0.05)*(self.filters["Wavelength_nm"]<900)
            
        myqe = self.qesmooth(self.filters["Wavelength_nm"][mask])
        return self.filters["Wavelength_nm"][mask]*10, self.filters[f][mask]*myqe/10000
        
    def get_avg_QE_filter(self, f):
        '''
        Returns the average value of the transmittance of the CCD due to filter response and QE.
        '''
        wl, resp = self.get_QE_filter( f)
        leff = self.lameff[f]
        dl = self.deltalam[f]

        return np.trapz((resp[np.abs(wl-leff)<dl/2]), wl[np.abs(wl-leff)<dl/2])  / dl
        
        
    
    def get_thoughput(self, log):
        
        cols = {"u":"purple", "r":"r", "g":"g", "i":"orange"}
        l = np.genfromtxt(log, dtype=None, delimiter=",", names=True)
        
        al_reflectivity = 0.86
        gain = 1.8
        
        avgth = {}   
        stdth = {}   
        
        plt.figure(figsize=(10,10))
        for i, f in enumerate(["u", "g", "r", "i"]):
            mask = l["filter"] ==f
            print (i%2,i/2)
            plt.subplot2grid((2,2), (i%2,i/2))
            
            qefilt = self.get_avg_QE_filter(f)
            
            obsflux = 10** (0.4 *( (2.5 * np.log10(l["exptime"][mask]) ) - l["inst"][mask])) * gain
            
            Ephot = (cnt.h).to(u.erg * u.s) * (cnt.c).to(u.angstrom/u.s) / (self.lameff[f]* u.angstrom)
            
            
            Fluxnu = (10 ** ( -0.4 * (l["std"][mask]+48.6 + self.get_Palomar_extinction(self.lameff[f], l["airmass"][mask]))) * u.erg * u.cm**-2 * u.s**-1 * u.Hz**-1).to(u.Jy)
            Flam = Fluxnu * cnt.c / (self.lameff[f]*u.angstrom)**2
            Flam =  Flam.to(u.erg /u.cm**2 / u.s/u.angstrom)
            
            predflux =  al_reflectivity * \
                qefilt * l["exptime"][mask] * u.s *\
                self.A * u.cm**2 *\
                (Flam / Ephot) * \
                self.lameff[f] * u.angstrom #* \
                #self.deltalam[f] * u.angstrom 
                #(((self.lameff[f]+self.deltalam[f]/2)* u.angstrom)**2 - ((self.lameff[f]-self.deltalam[f]/2)* u.angstrom)**2 )
                
            
            print "Flam",Flam
            print "AL_reflectivity %.2f QEFilt %.2f Area: %.2e cm**2 Energy photon %.2e erg"%(al_reflectivity, qefilt, self.A, Ephot.value)
            
            ratio = 1.* (obsflux/predflux).value
            avgth[f] = np.median(ratio)
            stdth[f] =  stats.funcs.median_absolute_deviation(ratio)

            mask2 = np.abs(ratio-np.median(ratio))< 5*stdth[f]
            
            plt.scatter(l["std"][mask][mask2], ratio[mask2], label="Predicted %s"%f, c=cols[f], alpha=0.5)
            #plt.ylim(0,1)
            
            print "Predicted", predflux[0:10],"observed",obsflux[0:10]
            plt.ylabel("Throughput (N$_{obs}$/N$_{pred}$)")
            plt.xlabel("PS1 / SDSS mag")
            
        
        plt.figure()
        for f in ["u", "g", "r", "i"]:
            plt.errorbar(self.lameff[f], avgth[f], yerr=stdth[f], fmt="o", color=cols[f], )
            
        plt.ylabel("Throughput (N$_{obs}$/N$_{pred}$)")
        plt.xlabel("Wavelength [$\\AA$]")
        plt.text(4000,0.05, "N($\\lambda$) = $\\frac{Alref QE(\\lambda) Filt(\\lambda) t_{exp}A 10^{-0.4(m_{AB} + 48.6 + EXT(\\lambda))} }{h\\lambda}$", fontsize=16)
        plt.show()
        
        return predflux
            
        
