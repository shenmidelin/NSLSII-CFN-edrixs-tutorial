# .bashrc

# Source global definitions
if [ -f /etc/bashrc ]; then
        . /etc/bashrc
fi

export PS1='$(whoami)@$(hostname):$(pwd)> '
export MY_SCRATCH=/hpcgpfs01/scratch/$(whoami)
alias squ="squeue -u $(whoami)"

module load intel/PSXE2018.u1
module load gcc/6.4.0
module load anaconda3/4.2.0

export LD_LIBRARY_PATH=/hpcgpfs01/work/workshop/edrixs_tutorial/lib/arpack_lib/lib:$LD_LIBRARY_PATH
export PYTHONPATH=/hpcgpfs01/work/workshop/edrixs_tutorial/lib/python3.5/site-packages:$PYTHONPATH
export PATH=/hpcgpfs01/work/workshop/edrixs_tutorial/bin:$PATH
