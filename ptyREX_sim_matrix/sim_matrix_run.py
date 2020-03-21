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
    """
    Running the matrix of sims.

    Parameters
    ----------
    xyz: str
    full path of the xyz coordination file
    sim_conditions: np.array
        array with the shape (3, n) with n the total number conditions under consideration
        The order of the three parameters:
        convergence semi-angle (rad), defocus (m), step_size (m)
    root_path: str
        full path of the directory where the sims are being saved
    script_path: str
        full path of the cluster script submit directory
    dose: int
        target sum intensity per diffraction pattern

    Returns
    -------
    """
    if not os.path.exists(root_path):
        os.makedirs(root_path)
    for i in np.arange(sim_conditions.shape[1]):
        sim_file = run_sim(script_path, root_path, xyz, sim_conditions[0,i], sim_conditions[1,i], sim_conditions[2,i])
        add_dose_noise(sim_file, dose)

    
if __name__ =='__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('xyz', help='path for the xyz file')
    parser.add_argument('sim_conditions', help='np.array with the sim conditions')
    parser.add_argument('root_path', help='path where the sims to be saved')
    parser.add_argument('script_path', help='path where the scripts live')
    parser.add_argument('dose', help='int, target sum dp')
    v_help = "Display all debug log messages"
    parser.add_argument("-v", "--verbose", help=v_help, action="store_true",
                        default=False)

    args = parser.parse_args()

    main(args.xyz, args.sim_conditions, args.root_path, args.script_path, args.dose)