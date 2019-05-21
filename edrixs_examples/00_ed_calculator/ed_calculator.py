#!/usr/bin/env python

"""
Warm up: use edrixs as a simple ED calculator
---------------------------------------------
Find all the eigenvalues and eigenvectors of a many-body Hamiltonian
with Coulomb interaction (and with spin-orbit coupling). Print out
the eigenvalues and their quantum numbers. 

You can change ``noccu`` from 0 to 6 to see what will happen.
"""

import numpy as np
import scipy
import edrixs

# set print options
np.set_printoptions(precision=6, linewidth=200, suppress=True)

# number of orbitals
norb= 6
# number of occupancy
noccu = 2
# Slater integrals
F0, F2 = 4.0, 1.0
# Coulomb interaction tensor
umat = edrixs.get_umat_slater('p', F0, F2)
# build Fock basis
basis = edrixs.get_fock_bin_by_N(norb, noccu)
# build Hamiltonian, four fermion terms
H = edrixs.build_opers(4, umat, basis)
# do ED
e1, v1 = scipy.linalg.eigh(H)
# add spin-orbit coupling (SOC)
soc = edrixs.atom_hsoc('p', 0.2)
# build Hamiltonian, two_fermion terms
H2 = H + edrixs.build_opers(2, soc, basis)
# do ED 
e2, v2 = scipy.linalg.eigh(H2)

orb_mom = edrixs.get_orb_momentum(1, True)
spin_mom = edrixs.get_spin_momentum(1)
tot_mom = orb_mom + spin_mom
opL, opS, opJ = edrixs.build_opers(2, [orb_mom, spin_mom, tot_mom], basis)
# L^2 = Lx^2 + Ly^2 +Lz^2, S^2 = Sx^2 + Sy^2 + Sz^2, J^2 = Jx^2 + Jy^2 + Jz^2
L2 = np.dot(opL[0], opL[0]) + np.dot(opL[1], opL[1]) + np.dot(opL[2], opL[2])
S2 = np.dot(opS[0], opS[0]) + np.dot(opS[1], opS[1]) + np.dot(opS[2], opS[2])
J2 = np.dot(opJ[0], opJ[0]) + np.dot(opJ[1], opJ[1]) + np.dot(opJ[2], opJ[2])

# print out results
print("Without SOC")
L2_val = edrixs.cb_op(L2, v1).diagonal().real
S2_val = edrixs.cb_op(S2, v1).diagonal().real
print("    #           eigenvalue     L(L+1)    S(S+1)")
for i, e in enumerate(e1):
    print("{:5d}{:20.6f}{:10.1f}{:10.1f}".format(i, e, L2_val[i], S2_val[i]))

print("With SOC")
print("    #           eigenvalue     J(J+1)")
J2_val = edrixs.cb_op(J2, v2).diagonal().real
for i, e in enumerate(e2):
    print("{:5d}{:20.6f}{:10.1f}".format(i, e, J2_val[i]))
