#!/usr/bin/env python

"""
Do ED of a system with 24 orbitals.

Usage: mpirun -np 4 python run_fedsolver.py [number_of_occupancy]
"""

import sys
import edrixs
from edrixs.fedrixs import ed_fsolver
from mpi4py import MPI

comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size()
fcomm = comm.py2f()

if len(sys.argv) > 1:
    noccu = int(sys.argv[1])
else:
    noccu = 12

if rank == 0:
    print("edrixs >>> build fock basis", noccu, "/", 24, flush=True)
    edrixs.write_fock_dec_by_N(24, noccu, "fock_i.in")
comm.Barrier()

ed_fsolver(fcomm, rank, size)
