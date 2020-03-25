#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
This is a script reads the pyprismatic sim parameters text file and creates
ptyREX json files to run reconstructions on DLS cluster.

"""

import numpy as np
import os
#from scipy import ndimage as ndi
import h5py
import argparse
import json
import collections
import sys
sys.path.append('/dls/science/groups/e02/Mohsen/code/Git_Repos/My_Repository/ptyREX_sim_matrix/')
from sim_utils import e_lambda

def parse_params_file(params_file, h5_file, drop_unneeded = True):
    '''
    Reads the parameters text file into a dict to be fed into ptypy / pycho recons
    '''
    exp_dict = {}
    with open(params_file) as f:    
        for line in f: 
            line = line.strip('-')
            exp_dict[line.strip().partition(':')[0]] = line.strip().partition(':')[-1]
    
    original_sim_file = exp_dict['output-file'].split('/')[-1]
    
    if original_sim_file != h5_file.split('/')[-1]:
        exp_dict['data'] = exp_dict.pop('output-file')
        exp_dict['data'] = h5_file
    else:
        exp_dict['data'] = exp_dict.pop('output-file')
        
    exp_dict['cell_dimension(m)'] = [1e-10*float(i) for i in exp_dict['cell-dimension'].split(' ')]
    exp_dict.pop('cell-dimension')
    exp_dict['tile_uc'] = [float(i) for i in exp_dict['tile-uc'].split(' ')]
    
    if 'skip' in h5_file:
        exp_dict['data_key'] = 'dataset'
    else:
        exp_dict['data_key'] = '4DSTEM_simulation/data/datacubes/hdose_noisy_data'
    exp_dict['xyz'] = exp_dict.pop('input-file')
    exp_dict['accel_voltage(eV)'] = int(exp_dict.pop('energy')) * 1000
    exp_dict['semi_angle(rad)'] = float(exp_dict.pop('probe-semiangle')) * 1e-3
    if 'skip' in h5_file:
        skip_ind = h5_file.index('skip')
        step_factor = int(h5_file[skip_ind + 4])
        exp_dict['step_size(m)'] = step_factor * float(exp_dict.pop('probe-step-x')) * 1e-10
    else:
        exp_dict['step_size(m)'] = float(exp_dict.pop('probe-step-x')) * 1e-10
    exp_dict['sim_pixel_size(m)'] = float(exp_dict.pop('pixel-size-x')) * 1e-10
    exp_dict['defocus(m)'] = float(exp_dict.pop('probe-defocus'))* 1e-10
    exp_dict['C3(m)'] = float(exp_dict.pop('C3'))* 1e-10
    exp_dict['C5(m)'] = float(exp_dict.pop('C5'))* 1e-10
    
    
    wavelength = e_lambda(exp_dict['accel_voltage(eV)'])
    
    exp_dict['wavelength'] = wavelength
    
    # getting rid of unneeded stuff:
    
    if drop_unneeded:
        exp_dict.pop('num-threads')
        exp_dict.pop('algorithm')
        exp_dict.pop('potential-bound')
        exp_dict.pop('num-FP')
        exp_dict.pop('slice-thickness')
        exp_dict.pop('num-slices')
        exp_dict.pop('zstart-slices')
        exp_dict.pop('alpha-max')
        exp_dict.pop('batch-size-cpu')
        exp_dict.pop('tile-uc')
        exp_dict.pop('detector-angle-step')
        exp_dict.pop('probe-xtilt')
        exp_dict.pop('probe-ytilt')
        exp_dict.pop('scan-window-x')
        exp_dict.pop('scan-window-y')
        exp_dict.pop('scan-window-xr')
        exp_dict.pop('scan-window-yr')
        exp_dict.pop('random-seed')
        exp_dict.pop('4D-amax')
        for k in ['thermal-effects', 'save-3D-output', 'save-4D-output', '4D-crop', 'save-DPC-CoM', 
                  'save-potential-slices','save-real-space-coords',  'occupancy', 'nyquist-sampling',
                  'probe-step-y', 'pixel-size-y']:
            exp_dict.pop(k)
    # using the sim parameters to calculate the bf disc rad
    det_pix_num = int((exp_dict['cell_dimension(m)'][0] * exp_dict['tile_uc'][0]) / (2 * exp_dict['sim_pixel_size(m)']))
    a_max = wavelength / (4 * exp_dict['sim_pixel_size(m)']) # alpha max
    print(a_max, wavelength)
    pix_per_rad = (det_pix_num / 2) / a_max
    exp_dict['pupil_rad(pixels)'] = pix_per_rad * exp_dict['semi_angle(rad)']
    
    exp_dict['detector_pixel_size(m)'] = 55e-6 # assuming the same as Medipix
    exp_dict['detector_distance(m)'] = (det_pix_num / 2) * exp_dict['detector_pixel_size(m)'] / a_max
    
    exp_dict['output_base'] = os.path.dirname(exp_dict['data'])
    
    exp_dict['rotation_angle(degrees)'] = 0
    
    return exp_dict
#%%


class NestedDefaultDict(collections.defaultdict):
    def __init__(self, *args, **kwargs):
        super(NestedDefaultDict, self).__init__(NestedDefaultDict, *args, **kwargs)

    def __repr__(self):
        return repr(dict(self))
    
def write_ptyrex_json(exp_dict):
    with h5py.File(exp_dict['data'], 'r') as f:
        data = f.get(exp_dict['data_key'])
        data_arr = np.array(data)
    
    scan_y = data_arr.shape[1]
    scan_x = data_arr.shape[0]
    
    N_x = data_arr.shape[2]
    N_y = data_arr.shape[3]
    
    # binning = np.floor(256 / N_x)
    
    # adj_px_size = exp_dict['detector_pixel_size(m)'] * binning

    

    params = NestedDefaultDict()
    
    params['process']['gpu_flag'] = 1
    params['process']['save_interval'] = 10
    params['process']['PIE']['iterations'] = 100
    params['process']['common']['source']['energy'] = [exp_dict['accel_voltage(eV)']]
    params['process']['common']['source']['radiation'] = 'electron'
    params['process']['common']['source']['flux'] = -1
    
    params['process']['common']['detector']['pix_pitch'] = list([exp_dict['detector_pixel_size(m)'], exp_dict['detector_pixel_size(m)']])
    params['process']['common']['detector']['distance'] = exp_dict['detector_distance(m)']
    params['process']['common']['detector']['bin'] = list([1, 1]) 
    params['process']['common']['detector']['min_max'] = list([0, 1000000])
    params['process']['common']['detector']['optic_axis']= list([N_x / 2, N_x/2])
    params['process']['common']['detector']['crop'] = list([N_x, N_y])
    params['process']['common']['detector']['orientation'] = '00'
    params['process']['common']['detector']['mask_flag'] = 0
    
    params['process']['common']['probe']['convergence'] = 2*exp_dict['semi_angle(rad)']
    params['process']['common']['probe']['distance'] = -1
    params['process']['common']['probe']['focal_dist'] = -1
    params['process']['common']['probe']['load_flag'] = 0
    params['process']['common']['probe']['diffuser'] = 0
    params['process']['common']['probe']['aperture_shape'] = 'circ'
    params['process']['common']['probe']['aperture_size'] = exp_dict['pupil_rad(pixels)']*exp_dict['detector_pixel_size(m)']

    params['process']['common']['object']['load_flag'] = 0
    
    params['process']['common']['scan']['rotation'] = exp_dict['rotation_angle(degrees)']
    params['process']['common']['scan']['fast_axis'] = 1
    params['process']['common']['scan']['orientation'] = '00'
    params['process']['common']['scan']['type'] = 'tv'
    params['process']['common']['scan']['load_flag'] = 0
    params['process']['common']['scan']['dR'] = list([exp_dict['step_size(m)'], exp_dict['step_size(m)']])
    params['process']['common']['scan']['N'] = list([scan_x, scan_y])
    
    params['experiment']['data']['data_path'] = exp_dict['data']
    params['experiment']['data']['dead_pixel_flag'] = 0 
    params['experiment']['data']['flat_field_flag'] = 0 
#    params['experiment']['data']['dead_pixel_path'] = exp_dict['mask']
#    params['experiment']['data']['flat_field_path'] = exp_dict['mask']
    params['experiment']['data']['load_flag'] = 1
    params['experiment']['data']['meta_type'] = 'hdf'
    params['experiment']['data']['key'] = exp_dict['data_key']
    
    params['experiment']['sample']['position'] = list([0, 0, 0])

    params['experiment']['detector']['position'] = list([0, 0, exp_dict['detector_distance(m)']])

    params['experiment']['optics']['lens']['alpha'] = 2*exp_dict['semi_angle(rad)']
    params['experiment']['optics']['lens']['defocus'] = list([exp_dict['defocus(m)'], exp_dict['defocus(m)']])
    params['experiment']['optics']['lens']['use'] = 1
    params['experiment']['optics']['diffuser']['use'] = 0
    params['experiment']['optics']['FZP']['use'] = 0
    params['experiment']['optics']['pinhole']['use'] = 0

    params['base_dir'] = exp_dict['output_base']
    params['process']['save_dir'] = exp_dict['output_base']
    params['process']['cores'] = 1
    
    json_file = os.path.join(exp_dict['output_base'], 'ptyREX_' + exp_dict['data'].split('/')[-1].split('.')[0] + '.json')
    exp_dict['ptyREX_json_file'] = json_file    
    with open(json_file, 'w+') as outfile:
        json.dump(params, outfile, indent = 4)


def main(params_file, h5_file):
    exp_dict = parse_params_file(params_file, h5_file)

    write_ptyrex_json(exp_dict)

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('params_txtFile', help='text file containing all the parameters used for simulation')
    parser.add_argument('h5_file', help ='full path of the h5 sim file')
    v_help = "Display all debug log messages"
    parser.add_argument("-v", "--verbose", help=v_help, action="store_true",
                        default=False)

    args = parser.parse_args()

    main(args.params_txtFile, args.h5_file)
