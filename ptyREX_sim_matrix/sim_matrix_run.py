#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 11:22:14 2020

@author: eha56862
"""
import os
import h5py
import pyprismatic as pr
import numpy as np
from shutil import copyfile
import argparse
import sys

#
#
#root_path = '/dls/e02/data/2020/cm26481-1/processing/pty_simulated_data_MD/sim_matrix_v3'
#
#if not os.path.exists(root_path):
#    os.mkdir(root_path)
#
#submit_path = '/dls/science/groups/e02/Mohsen/code/Git_Repos/My_Repository/create_sim_data_prismatic'
#
#
#coord_dict ={'/dls/science/groups/e02/Mohsen/code/Git_Repos/My_Repository/create_sim_data_prismatic/xyz_files/Graphene_SW_mod1.xyz':'graphene_small_hole'}
#             #'/dls/science/groups/e02/Mohsen/code/sim_4DSTEM/ptypy_pycho_sim_matrix/create_sim_data_prismatic/xyz_files/graphene_island_extended.xyz':'graphene_island',
#             #'/dls/science/groups/e02/Mohsen/code/sim_4DSTEM/ptypy_pycho_sim_matrix/create_sim_data_prismatic/xyz_files/graphene_island_doped_extended.xyz': 'graphene_island_doped'}
#
#convergence_dict = {10:'10mrad',
#                    15:'15mrad',
#                    25:'25mrad',
#                    32:'32mrad'}
#
#def_dict = {0:'zero_def',
#            50:'50A_def',
#            100:'100A_def',
#            150:'150A_def'}


def make_output_filename(root_path, xyz, conv_semiangle, def_val, step_size):
    '''
    makes an output filename -with path- that reflects the conditions used for sim
    Parameters
    ____________
    root_path: str
        full path of the directory where the sims are being saved
    xyz: str
        full path of the xyz coordination file used
    conv_semiangle: float
        probe convergence semi-angle (rad)
    def_val: float
        probe defocus (m)
    sep_size: float
        probe step size (m)
    Returns
    ______________
    output_file_path: str
        full path and h5 file sim file name
    '''
    output_name = os.path.basename(xyz).split('.')[0] + \
                '_' + str(conv_semiangle) + 'mrad_' + str(def_val) + \
                'm_def_' + str(step_size) + 'm_step_size'
    output_path = os.path.join(root_path, output_name)
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    
    output_file_path = os.path.join(output_path, output_name + '.h5')
    
    return output_file_path



def param_filename(root_path, xyz, conv_semiangle, def_val, step_size):
    '''
    makes an output parameters filename -with path- that reflects the conditions 
    used for sim
    Parameters
    ____________
    root_path: str
        full path of the directory where the sims are being saved
    xyz: str
        full path of the xyz coordination file used
    conv_semiangle: float
        probe convergence semi-angle (rad)
    def_val: float
        probe defocus (m)
    sep_size: float
        probe step size (m)
    Returns
    ______________
    output_file_path: str
        full path and txt file sim file name
    '''
    output_name = os.path.basename(xyz).split('.')[0] + \
                '_' + str(conv_semiangle) + 'mrad_' + str(def_val) + \
                'm_def_' + str(step_size) + 'm_step_size'
    output_path = os.path.join(root_path, output_name)
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    
    output_file_path = os.path.join(output_path, output_name + '.txt')
    
    return output_file_path



def copy_scratch_file(submit_path, root_path, xyz, conv_semiangle, def_val, step_size):
    '''
    Copies over the scratch file from the submission dir to corresponding sim dir
    Parameters
    ____________
    submit_path: str
        full path of the cluster script submit directory
    root_path: str
        full path of the directory where the sims are being saved
    xyz: str
        full path of the xyz coordination file used
    conv_semiangle: float
        probe convergence semi-angle (rad)
    def_val: float
        probe defocus (m)
    sep_size: float
        probe step size (m)
    Returns
    ______________

    
    '''
    scratch_file = os.path.join(submit_path, 'scratch_param.txt')
    copyfile(scratch_file, param_filename(root_path, xyz, conv_semiangle, def_val, step_size))
    
    

def get_cell_dims(xyz):
    '''
    returns the cell dimensions of an xyz atomic coordination file.
    Parameters
    ___________
    
    xyz: str
        full path of the xyz coordination file
    Returns
    ____________
    
    cell_dims: numpy array
        cell dimensions in (A)
    '''
    file = xyz
    with open(file, 'r') as f:
        for i, line in enumerate(f):
            if i == 1:
                data = np.asarray(line.split(), dtype = float)
    return data
    


def run_sim(submit_path, root_path, xyz, conv_semiangle, def_val, step_size):
    '''
    generates a meta parametre and runs a pyprismatic simulation
    
    Parameters
    _____________
    submit_path: str
        full path of the cluster script submit directory
    root_path: str
        full path of the directory where the sims are being saved
    xyz: str
        full path of the xyz coordination file
    conv_semiangle: float
        probe convergence semi angle in (rad)
    def_val: float
        defocus value in (m)
    ste_size: float
        probe step size in (m)
    Returns
    _____________
    
    '''
    
    meta = pr.Metadata(filenameAtoms = xyz)
    meta.algorithm = 'multislice'
    meta.filenameOutput = make_output_filename(root_path, xyz, conv_semiangle, def_val, step_size)
    sim_file_path = make_output_filename(root_path, xyz, conv_semiangle, def_val, step_size)
    # meta.writeParameters(param_filename(xyz, conv_semiangle, def_val))
    meta.numThreads = 12
    meta.realspacePixelSizeX = 0.11
    meta.realspacePixelSizeY = 0.11
    meta.potBound = 2
    meta.numFP = 8
    meta.sliceThickness = 8  # may change this 
    #meta.numSlices = 1
    #meta.zStart = 0
    meta.E0 = 80
    meta.alphaBeamMax = 45
    meta.batchSizeCPU = 1
    meta.probeStepX = step_size * 1e10
    meta.probeStepY = step_size * 1e10
    
    cell_dims = get_cell_dims(xyz)
    meta.cellDimX = cell_dims[0]
    meta.cellDimY = cell_dims[1]
    meta.cellDimZ = cell_dims[2]
    
#    if coord_dict[xyz] == 'graphene_bilayer':
#        meta.cellDimX = 16.7663
#        meta.cellDimY = 16.94
#        meta.cellDimZ = 3.395
#    elif coord_dict[xyz] == 'graphene_SW':
#        meta.cellDimX = 29.03
#        meta.cellDimY = 29.03
#        meta.cellDimZ = 1.1168
#    elif coord_dict[xyz] == 'graphene_hole':
#        meta.cellDimX = 81.5211
#        meta.cellDimY = 84.884
#        meta.cellDimZ = 8.000
    
    meta.tileX = 3
    meta.tileY = 3
    meta.tileZ = 1
    meta.probeDefocus = def_val * 1e10
    meta.C3 = 0
    meta.C5 = 0
    meta.probeSemiangle = conv_semiangle * 1e3
    meta.detectorAngleStep = 1
    meta.probeXtilt = 0
    meta.probeYtilt = 0
    meta.scanWindowXMin = 0.33
    meta.scanWindowXMax = 0.63
    meta.scanWindowYMin = 0.33
    meta.scanWindowYMax = 0.63
    #meta.scanWindowXMin_r = 0
    #meta.scanWindowXMax_r = 0
    #meta.scanWindowYMin_r = 0
    #meta.scanWindowYMax_r = 0
    meta.randomSeed = 25212
    meta.includeThermalEffects = 1
    meta.save2DOutput = 0
    meta.save3DOutput = 0
    meta.save4DOutput = 1
    meta.nyquistSampling = 0
    meta.saveDPC_CoM = 1
    meta.savePotentialSlices = 1
    meta.alsoDoCPUWork = 1
    meta.batchSizeGPU = 1
    meta.numGPUs = 2
    meta.numStreamsPerGPU = 3
    
    meta.go()
    
    copy_scratch_file(submit_path, root_path, xyz, conv_semiangle, def_val, step_size)
    
    return sim_file_path
    
def add_dose_noise(file_path, dose, add_noise = True):
    '''
    gets an h5 simulated 4DSTEM file and adds a dataset with the key:
        '4DSTEM_simulation/data/datacubes/hdose_noisy_data'
        in the original sim file, with dose multiplied to each frame and if passed True, 
        poisson noise.
    Parameters
    __________
    file_path: str
        full path and name of the sim h5 file
    dose: int
        target sum intensity per diffraction pattern
    add_noise: boolean
        if True it also adds posson noise to the diffraction patterns
    Returns
    ___________
    
    '''
    with h5py.File(file_path) as f:
        sh = f['4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube'].shape
        print('Dataset shape is %s' % str(sh))
        data = f.get('4DSTEM_simulation/data/datacubes/CBED_array_depth0000/datacube')
        data = np.array(data)
    
    
    if add_noise is False:
        data_highD = dose * data
    else:
        data_highD = dose * data
        data_highD = np.random.poisson(data_highD)
    f = h5py.File(file_path, 'a')
    f.create_dataset('4DSTEM_simulation/data/datacubes/hdose_noisy_data', data = data_highD, dtype='float32')
    f.close()
    
    return


#def save_skips(file_path):
#    with h5py.File(file_path) as f:
#        print(file_path)
#        sh = f['4DSTEM_simulation/data/datacubes/hdose_noisy_data'].shape
#        print('Dataset shape is %s' % str(sh))
#        data = f.get('4DSTEM_simulation/data/datacubes/hdose_noisy_data')
#        data = np.array(data)
#    skips = [2,3,4]
#    for skip in skips:
#        data_new = data[::skip, ::skip, :, :]
#        new_fn = os.path.basename(file_path).split('.')[0] + '_skip'+ str(skip)+ '.h5'
#        saving_path = os.path.join(os.path.dirname(file_path), new_fn)
#        hf = h5py.File(saving_path, 'w')
#        print('creating the h5 file for the stack')
#        hf.create_dataset('dataset', data=data_new, compression='gzip')
#    return
    

def main(xyz, sim_conditions, root_path, script_path, dose):
    
    
    for atom_model in list(coord_dict.keys()):
        for conv_semi in list(convergence_dict.keys()):
            for def_val in list(def_dict.keys()):
                sim_file = run_sim(atom_model, conv_semi, def_val)
                add_dose_noise(sim_file, 1e6)
#    
#    
#    for dirname, dirnames, filenames in os.walk(root_path):
#
#        for filename in filenames:
#            if filename.endswith('h5'):
#                save_skips(os.path.join(dirname, filename))

    
if __name__ =='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('beamline', help='Beamline name')
    parser.add_argument('year', help='Year')
    parser.add_argument('visit', help='Session visit code')
    parser.add_argument('folder', nargs= '?', help='Option to add folder')
    v_help = "Display all debug log messages"
    parser.add_argument("-v", "--verbose", help=v_help, action="store_true",
                        default=False)

    args = parser.parse_args()

    main(args.beamline, args.year, args.visit, args.folder)   