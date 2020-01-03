#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  2 11:22:14 2020

@author: eha56862
"""
import os
import pyprismatic as pr

root_path = '/dls/science/groups/e02/Mohsen/code/sim_4DSTEM/auto_sim/sim_matrix'


coord_dict ={'/dls/science/groups/e02/Mohsen/code/sim_4DSTEM/auto_sim/xyz_files/graphene_bilayers_4_7_tile.xyz':'graphene_bilayer',
             '/dls/science/groups/e02/Mohsen/code/sim_4DSTEM/auto_sim/xyz_files/Graphene_SW_cropped_square.xyz':'graphene_SW',
             '/dls/science/groups/e02/Mohsen/code/sim_4DSTEM/auto_sim/xyz_files/graphene_defect_new.xyz': 'graphene_hole'}

convergence_dict = {10:'10mrad',
                    15:'15mrad',
                    20:'20mrad',
                    25:'25mrad'}

def_dict = {0:'zero_def',
            100:'100A_def',
            200:'200A_def',
            300:'300A_def'}

def make_output_filename(xyz, conv_semiangle, def_val):
    '''
    makes an output filename -with path- that reflects the conditions used for sim
    '''
    output_name = coord_dict[xyz]+'_'+convergence_dict[conv_semiangle]+'_'+def_dict[def_val]
    output_path = os.path.join(root_path, output_name)
    if not os.path.exists(output_path):
        os.mkdir(output_path)
    
    output_file_path = os.path.join(output_path, output_name+'.h5')
    
    return output_file_path

def run_sim(xyz, conv_semiangle, def_val):
    '''
    generates a meta parametre and runs a pyprismatic simulation
    
    inputs:
        xyz: path of the xyz coordination file
        conv_semiangle: probe convergence semi angle in mrad
        def_val: defocus value in A
    '''
     
    meta = pr.Metadata(filenameAtoms = xyz)
    meta.algorithm = 'multislice'
    meta.filenameOutput = make_output_filename(xyz, conv_semiangle, def_val)
    meta.numThreads = 12
    meta.realspacePixelSizeX = 0.2
    meta.realspacePixelSizeY = 0.2
    meta.potBound = 2
    meta.numFP = 8
    meta.sliceThickness = 8  # may change this 
    #meta.numSlices = 1
    #meta.zStart = 0
    meta.E0 = 80
    meta.alphaBeamMax = 26
    meta.batchSizeCPU = 1
    meta.probeStepX = 0.2
    meta.probeStepY = 0.2
    
    if coord_dict[xyz] == 'graphene_bilayer':
        meta.cellDimX = 16.7663
        meta.cellDimY = 16.94
        meta.cellDimZ = 3.395
    elif coord_dict[xyz] == 'graphene_SW':
        meta.cellDimX = 29.03
        meta.cellDimY = 29.03
        meta.cellDimZ = 1.1168
    elif coord_dict[xyz] == 'graphene_hole':
        meta.cellDimX = 81.5211
        meta.cellDimY = 84.884
        meta.cellDimZ = 8.000
    
    meta.tileX = 1
    meta.tileY = 1
    meta.tileZ = 1
    meta.probeDefocus = def_val
    meta.C3 = 0
    meta.C5 = 0
    meta.probeSemiangle = conv_semiangle
    meta.detectorAngleStep = 1
    meta.probeXtilt = 0
    meta.probeYtilt = 0
    meta.scanWindowXMin = 0
    meta.scanWindowXMax = 0.99999
    meta.scanWindowYMin = 0
    meta.scanWindowYMax = 0.99999
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
    meta.savePotentialSlices = 0
    meta.alsoDoCPUWork = 1
    meta.batchSizeGPU = 1
    meta.numGPUs = 2
    meta.numStreamsPerGPU = 3
    
    meta.go()

def main():
    for atom_model in list(coord_dict.keys()):
        for conv_semi in list(convergence_dict.keys()):
            for def_val in list(def_dict.keys()):
                run_sim(atom_model, conv_semi, def_val)
    
if __name__ =='__main__':
    main()
    