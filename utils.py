#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:23:14 2020

@author: eha56862
"""
%matplotlib qt5
import numpy as np
import os
import hyperspy.api as hs
import h5py

def e_lambda(e_0):
    '''
    calculates the electron wavelength by taking the accelerating voltage in keV
    input:
        e_0 (accelerating voltage in keV)
    
    output:
        e_lambda (electron wavelength in A)
    '''
    import numpy as np
    
    m = 9.109383*10**-31 # electron rest mass in kg
    e = 1.602177*10**-19 # primary charge in C
    h = 6.62607*10**-34 # planck const in m**2*kg/s
    c =  299792458 # speed of light in m/s
    
    e_lambda = h/np.sqrt(2*m*e*e_0*10**3)/np.sqrt(1 + e*e_0*10**3/2/m/c**2)*10**10;
    
    return e_lambda

def sim_to_hs(sim_h5_file):
    '''
    reads simulated 4DSTEM file into hs object
    '''
    with h5py.File(sim_h5_file) as f:
        sh = f['4DSTEM_simulation/data/datacubes/hdose_noisy_data'].shape
        print('Dataset shape is %s' % str(sh))
        data = f.get('4DSTEM_simulation/data/datacubes/hdose_noisy_data')
        data = np.array(data)
        data_hs = hs.signals.Signal2D(data)
    return data_hs

def get_bf_disc(data_hs):
    cent_x = int(data_hs.axes_manager[3].size / 2)
    cent_y = int(data_hs.axes_manager[2].size / 2)
    circ_roi = hs.roi.CircleROI(cent_x,cent_y, 10)
    data_sum = data_hs.sum()
    data_sum.plot()
    imc = circ_roi.interactive(data_sum)

    return circ_roi
    


def get_haadf(data_hs, bf_rad):
    cent_x = int(data_hs.axes_manager[3].size / 2)
    cent_y = int(data_hs.axes_manager[2].size / 2)
    circ_roi = hs.roi.CircleROI(cent_x,cent_y, bf_rad+cent_x, bf_rad)
    data_T = data_hs.T
    data_T.plot()
    imc = circ_roi.interactive(data_T)
    imc.sum().plot()
    
    adf = imc.sum()
    return adf

def calc_camera_length(data_hs, bf_rad, angle, pixel_size):
    '''
    input:
        data_hs: 4D-STEM hypespy object
        bf_rad: bf radius in pixels
        angle: known angle in diff plane in mrad
        pixel_size: physcial size of detector pix in um
    returns
        CL: camera length in meters
    '''
    CL = bf_rad * pixel_size*1e-6 / (angle*1e-3)
    data_sum = data_hs.sum()
    cent_x = int(data_hs.axes_manager[3].size / 2)
    cent_y = int(data_hs.axes_manager[2].size / 2)
    circ_roi = hs.roi.CircleROI(cent_x,cent_y, bf_rad)
    data_sum.plot()
    circ_roi.interactive(data_sum)
    
    return CL

def make_mask(file_path, size):
    '''
    makes a transparent mask as h5 file with a given size
    '''
    mask = np.ones((size, size))
    with h5py.File(file_path, 'w') as f:
        dset = f.create_dataset('mask', data = mask)
    return
    
    
    
