# -*- coding: utf-8 -*-
"""
Created on Wed Jul 22 12:16:50 2015

@author: nadiablago
"""

from astropy.time import Time


def utc2mjd(times):
    t = Time(times, format='isot', scale='utc')
    return t.mjd
    
def mjd2utc(mjd):
    t = Time(mjd+2400000.5, format='jd', scale="utc")
    return t.iso
    
def jd2utc(jd):
    t = Time(jd, format='jd', scale="utc")
    return t.iso
    
def utc2jd(times):
    t = Time(times, format='iso', scale='utc')
    return t.jd