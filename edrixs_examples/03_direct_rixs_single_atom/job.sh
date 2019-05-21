#!/bin/bash

#SBATCH --job-name=Ni_rixs
#SBATCH -p long
#SBATCH -t 00:30:00
#SBATCH -A workshop
#SBATCH -N 1
#SBATCH -n 1
#SBATCH --qos normal

module load intel/PSXE2018.u1
module load gcc/6.4.0
module load anaconda3/4.2.0

export LD_LIBRARY_PATH=/hpcgpfs01/work/workshop/edrixs_tutorial/lib/arpack_lib/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/hpcgpfs01/work/workshop/edrixs_tutorial/lib/python3.5/site-packages:$PYTHONPATH

srun -N 1 -n 1 python run_rixs_pysolver.py > log.txt
