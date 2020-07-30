#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jun 22 23:01:50 2020

@author: eha56862
"""

excel_path = '/dls/e02/data/2020/cm26481-1/processing/Merlin/20200130_80kV_graphene_600C_pty/20200130_Ptychography_speadsheet.xlsx'
session_path = '/dls/e02/data/2020/cm26481-1/processing/Merlin/20200130_80kV_graphene_600C_pty/'
header = 7
cols = "A:R" 

#%%
import pandas as pd
import os
import h5py
import pyxem as pxm
import numpy as np
import sys
sys.path.append('/dls/science/groups/e02/Mohsen/code/Git_Repos/Merlin-Medipix/')
import epsic_tools.api as epsic

#%%
def parse_excel(excel_path):
    data = pd.read_excel(excel_path, header=header, usecols=cols)
    # remoce rows where all elements are nans
    data.dropna(how='all', inplace=True) 
    
    list_of_dicts = []
    for i in range(len(data)):
        entry = data.iloc[i]
        entry_dict = entry.to_dict()
        list_of_dicts.append(entry_dict)
    return list_of_dicts

#%%
def _create_metadata_dict(session_path, entry_dict):
    metadata_dict = {}
    
    metadata_dict['dataset'] = entry_dict['Folder']
    metadata_dict['accelerating voltage (keV)'] = entry_dict.pop('keV')
    metadata_dict['A2 value (keV)'] = entry_dict.pop('A2 (kV)')
    metadata_dict['spot size'] = entry_dict.pop('Spot')
    metadata_dict['CL aperture'] = entry_dict.pop('CLA')
    metadata_dict['convergence simi-angle'] = {}
    metadata_dict['convergence simi-angle']['nominal (mrad)'] = entry_dict.pop('alpha')
    metadata_dict['convergence simi-angle']['calibrated (mrad)'] = []
    metadata_dict['convergence simi-angle']['calibration comments'] = []
    metadata_dict['camera length'] = {}
    metadata_dict['camera length']['nominal (cm)'] = entry_dict.pop('CL')
    metadata_dict['camera length']['calibrated (cm)'] = []
    metadata_dict['camera length']['calibration comments'] = []
    metadata_dict['defocus'] = {}
    metadata_dict['defocus']['nominal (nm)'] = entry_dict.pop('nominal defocus (nm)')
    metadata_dict['defocus']['calibrated (nm)'] = []
    metadata_dict['defocus']['calibration comments'] = []
    metadata_dict['MAG'] = entry_dict.pop('Scan Mag')
    metadata_dict['FOV (A)'] = entry_dict.pop('Field of view (Ang)')
    metadata_dict['recorded probe positions'] = entry_dict.pop('probe positions recorded')
    metadata_dict['counter depth'] = entry_dict.pop('counter depth')
    metadata_dict['saved data bit depth'] = entry_dict.pop('Saved data bit depth')
    metadata_dict['step size (A)'] = entry_dict.pop('step size (Ang)')
    metadata_dict['frame time (ms)'] = entry_dict.pop('frame time (ms)')
    
    metadata_dict['Comments'] = entry_dict.pop('Notes')
    
    return metadata_dict
#%% 
#test = '/dls/e02/data/2020/cm26481-1/processing/Merlin/20200130_80kV_graphene_600C_pty/20200131 115200/80kV_600C_CLA_40um_CL_8cm_8C_20Mx_A2_4p71_df0nm_scan_array_255by255_diff_plane_515by515_.hdf5'

def _get_hdf_details(hdf_file):
    data = pxm.load(hdf_file, lazy=True)
    return data.data.shape


def _get_hdf_list(path):
    hdf_list = []
    for file in os.listdir(path):
        if file.endswith('.hdf5'):
            hdf_list.append(os.path.join(path,file))
    return hdf_list

        
#%%
list_of_dicts = parse_excel(excel_path)
#%%
for entry in list_of_dicts:
    results = []
    if entry['Folder'] in os.listdir(session_path):
        metadata_dict = _create_metadata_dict(session_path, entry)
        results.append(metadata_dict)
        hdf5_files = _get_hdf_list(os.path.join(session_path, entry['Folder']))
        for file in hdf5_files:
            a, b = os.path.splitext(file)
            meta_file_name = a + '_metadata.hdf5'
            epsic.ptycho_utils.save_dict_to_hdf5(metadata_dict, meta_file_name)
#            save_metadata_to_hdf5(metadata_dict, meta_file_name)
    else:
        print(entry['Folder'] , ' not among the datasets - no metadata added.')

#%%
#example dict I am writing:
metadata_ex = {'dataset': '20200130 150828',
  'accelerating voltage (keV)': 80.0,
  'A2 value (keV)': 4.71,
  'spot size': '8C',
  'CL aperture': '40um',
  'convergence simi-angle': {'nominal (mrad)': 31.74,
   'calibrated (mrad)': [],
   'calibration comments': []},
  'camera length': {'nominal (cm)': 33.616,
   'calibrated (cm)': [],
   'calibration comments': []},
  'defocus': {'nominal (nm)': 0.0,
   'calibrated (nm)': [],
   'calibration comments': []},
  'MAG': '20Mx',
  'FOV (A)': 9.52e-09,
  'recorded probe positions': 256.0,
  'counter depth': 6.0,
  'saved data bit depth': 8.0,
  'step size (A)': 3.71875e-11,
  'frame time (ms)': 0.7,
  'Comments': 'CL - 20cm on microscope,  sampling of probe too low as  bfd takes up more than half of detector '}
