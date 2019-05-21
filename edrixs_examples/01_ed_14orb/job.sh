#!/bin/bash

#SBATCH --job-name=ed_14orb
#SBATCH -p long
#SBATCH -t 00:30:00
#SBATCH -A workshop
#SBATCH -N 1
#SBATCH -n 8
#SBATCH --qos normal

module load intel/PSXE2018.u1
module load gcc/6.4.0
module load anaconda3/4.2.0

export LD_LIBRARY_PATH=/hpcgpfs01/work/workshop/edrixs_tutorial/lib/arpack_lib/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/hpcgpfs01/work/workshop/edrixs_tutorial/lib/python3.5/site-packages:$PYTHONPATH
export PATH=/hpcgpfs01/work/workshop/edrixs_tutorial/bin:$PATH

# option 1: get input files and also do ED with pure Python solver
srun -N 1 -n 1 python run_ed_pysolver.py > log1.txt

# option 2: call standalone executables of Fortran ED solver: ed.x 
srun -N 1 -n 8 ed.x > log2.txt
cp eigvals.dat eigvals_edx.dat

# option 3: use Python API: ed_fsolver to call Fortran ED solver
srun -N 1 -n 8 python run_ed_fsolver.py > log3.txt
