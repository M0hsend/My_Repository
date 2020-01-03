#!/bin/bash

#$ -P tomography
#$ -N pyprismatic_epsic
#$ -l gpu_arch=Pascal
#$ -o stdout_mc2.log -e stderr_mc2.log
#$ -l exclusive
#$ -q high.q
#$ -l gpu=2
#$ -cwd

module load python/epsic3.7


python sim_matrix.py
