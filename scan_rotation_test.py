#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb  7 12:10:14 2020

This script gets the path of json file and generates a series of json files with 
inceremental rotation angles in scan to check the effect on pty recon

@author: eha56862
"""

import json
import os
import scandir
import argparse


json_file = '/dls/e02/data/2020/cm26481-1/processing/pty_simulated_data_MD/sim_matrix_v2/graphene_small_hole_20mrad_200A_def/rotation_test/ptyREX_graphene_small_hole_20mrad_200A_def.json'

def write_json_series(json_file, rot_inc, steps_num):
    with open(json_file) as r:
        params = json.load(r)
    
    rot_start = params['process']['common']['scan']['rotation']
    rots = []
    for i in range(steps_num):
        rots.append(rot_start + rot_inc * i)
    for new_rot in rots:
        params['process']['common']['scan']['rotation'] = new_rot
        save_dir = os.path.join(os.path.dirname(json_file), os.path.splitext(os.path.basename(json_file))[0]+ '_'+ str(new_rot)+'_deg')
        if not os.path.exists(save_dir):
            os.mkdir(save_dir)
        params['base_dir'] = save_dir
        params['process']['save_dir'] = save_dir
        
        json_file_new = os.path.join(save_dir, os.path.splitext(os.path.basename(json_file))[0]+ '_'+ str(new_rot)+'_deg'+'.json')
        
        with open(json_file_new, 'w+') as outfile:
            json.dump(params, outfile, indent = 4)
    return

def get_json_list(sim_matrix_path):
    '''
    checks for the folders that have json files
    '''
    json_dirs = []
    it =  scandir.scandir(sim_matrix_path)
    for entry in it:
        if entry.is_dir():
            it2 = scandir.scandir(entry.path)
            for entry2 in it2:
                if entry2.is_file():
                    if entry2.name.endswith('.json'):
                        json_dirs.append(os.path.join(entry.path, entry.name))
    return json_dirs


def main(json_path, rot_inc, steps_num):
    script_path = '/dls/science/groups/e02/Mohsen/code/sim_4DSTEM/ptypy_pycho_sim_matrix/'
    write_json_series(json_path, rot_inc, steps_num)
    base_dir = os.path.dirname(json_path)
    json_list = get_json_list(base_dir)
    for json_file in json_list:
        json_noext = os.path.splitext(json_file)[0]
        output_folder = os.path.dirname(json_file)
        os.system('\n cd '+ output_folder + '\n module load global/cluster \n qsub '+ script_path + 'ptyREX_batch_submit.sh '+ output_folder + ' ' + json_noext + ' 07022020')
    

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('json_path', help='path to starting json file')
    parser.add_argument('rot_inc', type = int, help='step applied')
    parser.add_argument('step_num', type = int, help='number of rotation steps')
    v_help = "Display all debug log messages"
    parser.add_argument("-v", "--verbose", help=v_help, action="store_true",
                        default=False)

    args = parser.parse_args()

    main(args.json_path, args.rot_inc, args.step_num)