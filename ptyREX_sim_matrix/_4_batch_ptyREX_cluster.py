#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 22 11:19:23 2020

@author: eha56862
This script gets the folder path with a sim matrix and does the following:
     - checks which folders have seen no action at all, i.e. only contain two files - .h5 data file and .txt param file
     - runs the 4_data_input_prep_ptyREX.py on the folder above to prepare them for recons
     - submits ptyREX jobs from these folders 
"""

import os
import sys
import numpy as np
import argparse


def get_raw_dir_list(sim_matrix_path, get_all = False):
    '''
    checks for the folders with only two files and identify them as raw
    get_all set to True returns all the folders 
    '''
    raw_dirs = []
    it =  os.scandir(sim_matrix_path)
    if get_all:
        for entry in it:
            if entry.is_dir():
                raw_dirs.append(entry.path)
    else:
        for entry in it:
            if entry.is_dir():
     #           if 'recons' not in os.listdir(os.path.dirname(entry.path)): 
                if len(os.listdir(entry.path)) == 2:
                    raw_dirs.append(entry.path)
    return raw_dirs


def get_ptyREX_ready(sim_matrix_path):
    '''
    checks for the folders that have ptyREX json file 
    '''
    ptyREX_dirs = []
    it =  os.scandir(sim_matrix_path)
    for entry in it:
        if entry.is_dir():
            it2 = os.scandir(entry.path)
            for entry2 in it2:
                if entry2.is_file():
                    if entry2.name.startswith('ptyREX_'):
                        ptyREX_dirs.append(entry.path)
    return ptyREX_dirs
    

def main(sim_matrix_path):
    script_path = '/dls/science/groups/e02/Mohsen/code/Git_Repos/My_Repository/ptyREX_sim_matrix'
#    dirs_to_prep = get_raw_dir_list(sim_matrix_path, get_all = True)
#    for path in dirs_to_prep:
#        it = os.scandir(path)
#        for entry in it:
#            if entry.is_file():
#                if entry.name.endswith('txt'):
#                    entry_path = '/'+os.path.join(*entry.path.split('/')[:-1])
#                    for item in os.listdir(entry_path):
#                        if item.endswith('h5'):
#                            print(os.path.join(entry_path, item))
#                            os.system('\n python '+ script_path+'/_4_data_input_prep_ptyREX.py '+ entry.path + ' ' + os.path.join(entry_path, item))
#                            os.system('\n module load global/cluster \n qsub '+ script_path + '/_4_data_prep_submit.sh '+ entry.path+ ' '+ item)

    dirs_to_run_ptyREX = get_ptyREX_ready(sim_matrix_path) 
    for path in dirs_to_run_ptyREX:
        it = os.scandir(path)
        for entry in it:
            if entry.is_file():
                if entry.name.startswith('ptyREX'):
                    output_folder = os.path.dirname(entry.path)
                    json_file = os.path.splitext(entry.name)[0]
                    os.system('\n cd '+ output_folder + '\n module load global/cluster \n qsub '+ script_path + '/ptyREX_batch_submit.sh '+ output_folder + ' ' + json_file + ' 25032020')
        
        
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('sim_matrix_path', help='path containing all the simulated sets')
    v_help = "Display all debug log messages"
    parser.add_argument("-v", "--verbose", help=v_help, action="store_true",
                        default=False)

    args = parser.parse_args()

    main(args.sim_matrix_path)
