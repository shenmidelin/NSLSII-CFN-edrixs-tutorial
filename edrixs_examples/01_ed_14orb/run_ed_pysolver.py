#!/usr/bin/env python

import numpy as np
import scipy
import edrixs

if __name__ == "__main__":
    """
    Do exact diagonalization of a f electronic system with 14 orbitals.

    You can change noccu to play with it.
    """
    # number of orbitals
    norbs = 14

    # number of occupancy, change this number and re-run it to see how it will change.
    noccu = 7

    # Slater integrals: FX, x=0, 2, 4, 2*3=6
    F2_f, F4_f, F6_f = 9.711 * 0.77, 6.364 * 0.77, 4.677 * 0.77
    # set the average Coulomb interaction energy Uf_ave to zero
    Uf_av = 0.0
    F0_f = Uf_av + edrixs.get_F0('f', F2_f, F4_f, F6_f)
    # get rank-4 Coulomb U-tensor
    params = [F0_f, F2_f, F4_f, F6_f]
    umat_i = edrixs.get_umat_slater('f', *params)
    # SOC strength
    zeta_f_i = 0.261 * 0.9
    hsoc_i = edrixs.atom_hsoc('f', zeta_f_i)


    # prepare files for ed.x
    # write control parameters to file 
    edrixs.write_config(ed_solver=2, num_val_orbs=norbs, neval=100, ncv=200, nvector=1, idump=True)
    # write nonzeros terms of two-fermion terms hsoc_i to file 'hopping_i.in', read by ed.x
    edrixs.write_emat(hsoc_i, 'hopping_i.in', 1E-10)
    # write nonzeros terms of four-fermion terms umat to file 'coulomb_i.in', read by ed.x
    edrixs.write_umat(umat_i, 'coulomb_i.in', 1E-10)
    # write fock basis of decimal format to file 'fock_i.in', read by ed.x
    edrixs.write_fock_dec_by_N(norbs, noccu, "fock_i.in")


    # we also use pure Python ED solver to get the eigenvalues
    print("edrixs >>> building Fock basis ...")
    basis_i = edrixs.get_fock_bin_by_N(norbs, noccu)
    print("edrixs >>> Done!")

    print("edrixs >>> building Hamiltonian ...")
    H = edrixs.build_opers(2, hsoc_i, basis_i)
    H += edrixs.build_opers(4, umat_i, basis_i)
    print("edrixs >>> Done!")

    print("edrixs >>> diagonalize Hamiltonian ...")
    eval_i, evec_i = scipy.linalg.eigh(H)
    edrixs.write_tensor(eval_i, "eval_i.dat", fmt_float='{:20.10f}')
    print("edrixs >>> Done!")
