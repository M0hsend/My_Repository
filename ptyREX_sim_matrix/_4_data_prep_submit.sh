#!/bin/bash
#$ -l h_rt=30:00:00
#$ -q high.q
#$ -l redhat_release=rhel7
#$ -l m_mem_free=100G
#$ -o /dls/science/groups/e02/Mohsen/code/Git_Repos/My_Repository/ptyREX_sim_matrix/logs/logs_prep
#$ -e /dls/science/groups/e02/Mohsen/code/Git_Repos/My_Repository/ptyREX_sim_matrix/logs/logs_prep

python _4_data_input_prep_ptyREX.py $1 $2
