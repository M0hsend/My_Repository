#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Feb 11 14:23:14 2020

@author: eha56862
"""
# %matplotlib qt5
import numpy as np
import hyperspy.api as hs
import h5py
from scipy import constants as pc


def e_lambda(e_0):
    """
    relativistic electron wavelength

    :param e_0: int
        accelerating voltage in volts
    :return:
    e_lambda: float
        wavelength in meters
    """
    import numpy as np
    
    # m = 9.109383*10**-31 # electron rest mass in kg
    # e = 1.602177*10**-19 # primary charge in C
    # h = 6.62607*10**-34 # planck const in m**2*kg/s
    # c =  299792458 # speed of light in m/s
    
    e_lambda = (pc.h * pc.c) / np.sqrt((pc.e * e_0)**2  + 2 * pc.e * e_0 * pc.m_e * pc.c**2)
    
    return e_lambda


def sim_to_hs(sim_h5_file, h5_key = 'hdose_noisy_data'):
    '''
    reads simulated 4DSTEM file into hs object
    '''
    if h5_key == 'hdose_noisy_data':
        with h5py.File(sim_h5_file, 'r') as f:
            sh = f['4DSTEM_simulation/data/datacubes/hdose_noisy_data'].shape
            print('Dataset shape is %s' % str(sh))
            data = f.get('4DSTEM_simulation/data/datacubes/hdose_noisy_data')
            data = np.array(data)
            data_hs = hs.signals.Signal2D(data)
    elif h5_key == 'skip':
        with h5py.File(sim_h5_file, 'r') as f:
            sh = f['dataset'].shape
            print('Dataset shape is %s' % str(sh))
            data = f.get('dataset')
            data = np.array(data)
            data_hs = hs.signals.Signal2D(data)
        
    else:
        with h5py.File(sim_h5_file, 'r') as f:
            sh = f['4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube'].shape
            print('Dataset shape is %s' % str(sh))
            data = f.get('4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube')
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
        angle: known angle in diff plane in rad
        pixel_size: physcial size of detector pix in m
    returns
        CL: camera length in meters
    '''
    CL = bf_rad * pixel_size / angle
    data_sum = data_hs.sum()
    cent_x = int(data_hs.axes_manager[3].size / 2)
    cent_y = int(data_hs.axes_manager[2].size / 2)
    circ_roi = hs.roi.CircleROI(cent_x,cent_y, bf_rad)
    data_sum.plot()
    circ_roi.interactive(data_sum)
    
    return CL

    
def get_disc_overlap(rad, dist):
    """
    to calculate disc over lap in image or diff plane
    
    rad: float
        radius of probe or probe semi-angle
    dist: float
        step distance or the angle to the reflection of interest
    returns
    
    Percentage_overlap: float
        
    """
    #r =3.2e-10#0.0155 #1e-10 #probe radius or convergence semi angle
    #d = 4*4.5e-11#0.0196#4*3.72e-11# probe spacing or reflection in mrad (atomic_spacing_rad)
    x_pos  = 0.5 * dist # x coordinate of circle intersection
    y_pos = np.sqrt(rad**2 - x_pos**2) # y coordinate of circle intersection
    theta = 2*np.arctan(y_pos / x_pos) # angle subtended by overlap
    A_overlap =  (rad**2 * theta) - (2*x_pos*y_pos) # area of overlap
    A_probe = np.pi * rad**2
    Percentage_overlap =100 *  A_overlap / A_probe
    print('x, y : ', x_pos, y_pos)
    print('theta : ', theta )
    print('overlap area : ', A_overlap)
    print('probe area : ', A_probe)
    print('overlap  % : ', Percentage_overlap)
    return Percentage_overlap

def get_overlap(probe_rad, step_size):
    """
    probe_rad: float
        probe radius in A
    step_size: float
        scan step size
    Returns
    probe_overlap: float
        percentage probe overlap
    """
    probe_overlap = 1 - step_size / (2 * probe_rad)
    
    return 100 * probe_overlap


def get_step_size(probe_rad, target_overlap):
    """
    knwoing the probe radius and the target overlap percentage this function returns
    the suitable step size.
    Parameters
    ___________
    probe_rad: float
        probe radius in A
    target_overlap: float
        overlap fraction
    Returns
    _________
    step_size: float
        the step size needed to get the target overlap
    """
    step_size = (1 - target_overlap) * (2 * probe_rad)
    
    return step_size