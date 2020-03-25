#!/bin/bash

#$ -P tomography
#$ -N pyprismatic_epsic
#$ -l gpu_arch=Pascal
#$ -o /dls/science/groups/e02/Mohsen/code/Git_Repos/My_Repository/ptyREX_sim_matrix/logs/logs_pyprismatic
#$ -e /dls/science/groups/e02/Mohsen/code/Git_Repos/My_Repository/ptyREX_sim_matrix/logs/logs_pyprismatic
#$ -l exclusive
#$ -l m_mem_free=150G
#$ -q high.q
#$ -l gpu=2
#$ -cwd

module load python/epsic3.7


python 2_sim_matrix_run.py
