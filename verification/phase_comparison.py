#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 10 16:33:53 2020

Gets the potential output from the pyprismatic sim and outputs the expected phase shift

@author: eha56862
"""
import numpy as np
import h5py

def get_potential(sim_file_path):
    '''
    gets the pyprismatic h5 file and outputs the potential - in V.Angstroms
    '''
    with h5py.File(sim_file_path, 'r') as f:
        pots = f['4DSTEM_simulation']['data']['realslices']['ppotential']['realslice'][:]
        pots = np.squeeze(pots)
    return pots

def e_lambda(e_0):
    """
    calculates the electron wavelength by taking the accelerating voltage in keV
    input:
        e_0 (accelerating voltage in keV)
    
    output:
        e_lambda (electron wavelength in A)
    """
    import numpy as np
    
    m = 9.109383*10**-31 # electron rest mass in kg
    e = 1.602177*10**-19 # primary charge in C
    h = 6.62607*10**-34 # planck const in m**2*kg/s
    c =  299792458 # speed of light in m/s
    
    e_lambda = h/np.sqrt(2*m*e*e_0*10**3)/np.sqrt(1 + e*e_0*10**3/2/m/c**2)*10**10;
    
    return e_lambda

def _sigma(e_0):
    """
    From Pete Nellist MATLAB code
    return the interaction parameter sigma in radians/(Volt-Angstroms)
    ref: Physics Vade Mecum, 2nd edit, edit. H. L. Anderson
    The American Institute of Physics, New York) 1989
     page 4.

    :param e_0: accelerating voltage in keV
    :return: sigma - the interaction parameter sigma in radians/(Volt-Angstroms)
    """
    emass = 510.99906 #electron rest mass in keV
    l = e_lambda(e_0)
    x = (emass + e_0) / (2.0 * emass + e_0)
    s = 2.0 * np.pi * x / (l * e_0)
    s = s / 1000 # in radians / (V.A)
    return s

