#!/bin/bash

#SBATCH --job-name=all_job
#SBATCH -p long
#SBATCH -t 00:30:00
#SBATCH -A workshop
#SBATCH -N 1
#SBATCH -n 4
#SBATCH --qos normal

module load intel/PSXE2018.u1
module load gcc/6.4.0
module load anaconda3/4.2.0

export LD_LIBRARY_PATH=/hpcgpfs01/work/workshop/edrixs_tutorial/lib/arpack_lib/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/hpcgpfs01/work/workshop/edrixs_tutorial/lib/python3.5/site-packages:$PYTHONPATH
export PATH=/hpcgpfs01/work/workshop/edrixs_tutorial/bin:$PATH

cd ed
srun -N 1 -n 4 ed.x > log.txt
cp eigvec.* ../xas
cp eigvec.* ../rixs_pp
cp eigvec.* ../rixs_ps

cd ../xas
srun -N 1 -n 4 xas.x > log.txt

cd ../rixs_pp
srun -N 1 -n 4 rixs.x > log.txt

cd ../rixs_ps
srun -N 1 -n 4 rixs.x > log.txt

# plot xas and rixs spectra
cd ..
get_spectrum.py -N 1000 -ommin -20 -ommax 10 -G 2.5 -T 10 -off 11219 -f xas.dat xas/xas_poles*
get_spectrum.py -N 1000 -ommin -0.2 -ommax 1.5 -G 0.04 -T 10 -f rixs.dat rixs_pp/rixs_poles* rixs_ps/rixs_poles*
