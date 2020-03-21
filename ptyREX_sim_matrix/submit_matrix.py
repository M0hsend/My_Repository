import numpy as np
import os

sim_conds = np.load(sim_conditions)
xyz = '/dls/science/groups/e02/Mohsen/code/Git_Repos/My_Repository/ptyREX_sim_matrix/xyz_files/Graphene_defect.xyz'
root_path = '/dls/e02/data/2020/cm26481-1/processing/pty_simulated_data_MD/sim_matrix_ptyREX/'
script_path = '/dls/science/groups/e02/Mohsen/code/Git_Repos/My_Repository/ptyREX_sim_matrix/'
dose = 1000

os.system('\n cd ' + script_path+'\n')
os.system('module load global/cluster \n')
os.system('\n qsub prismatic_matrix.sc ' + xyz + ' ' + sim_conds + ' ' + root_path + ' ' + script_path + ' ' + dose)