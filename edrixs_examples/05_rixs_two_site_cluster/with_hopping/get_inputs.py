#!/usr/bin/env python


import os
import shutil
import numpy as np
import edrixs


def build_dirs():
    for path in ["ed", "xas", "rixs_pp", "rixs_ps"]:
        if os.path.isdir(path):
            shutil.rmtree(path)
        os.mkdir(path)
    shutil.copy('job_ed.sh',     'ed/job_ed.sh')
    shutil.copy('job_xas.sh',   'xas/job_xas.sh')
    shutil.copy('job_rixs.sh', 'rixs_pp/job_rixs.sh')
    shutil.copy('job_rixs.sh', 'rixs_ps/job_rixs.sh')


def get_hopping_coulomb(locaxis):
    # Number of orbitals for each site
    ndorb, nporb = 6, 4
    # Number of sites
    nsite = 2
    # Total number of orbitals
    ntot = nsite * (ndorb + nporb)
    # orbital orders:
    # 0-5:    1st-site-t2g
    # 6-11:   2nd-site-t2g
    # 12-15:  1st-site-2p
    # 16-19:  2nd-site-2p

    # On-site Coulomb interaction tensor
    U, J = 2.0, 0.3  
    Ud, JH = edrixs.UJ_to_UdJH(U, J)
    F0_dd, F2_dd, F4_dd = edrixs.UdJH_to_F0F2F4(Ud, JH)   # k=0, 2, 2*l

    G1_dp, G3_dp = 0.957, 0.569  # k=|2-1|, |2+1|
    F0_dp, F2_dp = edrixs.get_F0('dp', G1_dp, G3_dp), 1.107  # k=0, min(2*2, 2*1)

    # just one site t2g-subspace
    umat_tmp_i = edrixs.get_umat_slater('t2g', F0_dd, F2_dd, F4_dd)

    params = [
        F0_dd, F2_dd, F4_dd,    # FX for dd
        F0_dp, F2_dp,           # FX for dp
        G1_dp, G3_dp,           # GX for dp
        0, 0                    # FX for pp
    ]
    # just one site
    umat_tmp_n = edrixs.get_umat_slater('t2gp32', *params)  # 2p_3/2 -> t2g

    # static core-hole potential
    static_v = 2.0
    for i in range(0, ndorb):
        for j in range(ndorb, ndorb + nporb):
            umat_tmp_n[i, j, j, i] += static_v

    # two sites as a whole
    umat_i = np.zeros((ntot, ntot, ntot, ntot), dtype=np.complex)
    umat_n = np.zeros((ntot, ntot, ntot, ntot), dtype=np.complex)

    umat_i[0:ndorb, 0:ndorb, 0:ndorb, 0:ndorb] = umat_tmp_i  # 1st site 5d-valence
    umat_i[ndorb:2*ndorb, ndorb:2*ndorb, ndorb:2*ndorb, ndorb:2*ndorb] = umat_tmp_i  # 2nd site 5d-valence

    indx = np.array([[ 0,  1,  2,  3,  4,  5,  # orbital indices for 1st site 5d-t2g
                      12, 13, 14, 15],         # orbital indices for 1st site 2p-core
                     [ 6,  7,  8,  9, 10, 11,  # orbital indices for 2nd site 5d-t2g
                      16, 17, 18, 19]          # orbital indices for 2nd site 2p-core
                    ])
    # copy umat_tmp_n (one site) to umat_n (two sites)
    ndp = ndorb + nporb
    for m in range(nsite):
        for i in range(ndp):
            for j in range(ndp):
                for k in range(ndp):
                    for l in range(ndp):
                        umat_n[indx[m, i], indx[m, j], indx[m, k], indx[m, l]] += umat_tmp_n[i, j, k, l]

    # two fermion terms, SOC, crystal field, and hopping between the two sites
    emat_i = np.zeros((ntot, ntot), dtype=np.complex)
    emat_n = np.zeros((ntot, ntot), dtype=np.complex)

    # SOC
    zeta_d_i = 0.35
    soc_d = edrixs.atom_hsoc('t2g', zeta_d_i)

    emat_i[0:ndorb, 0:ndorb] += soc_d
    emat_i[ndorb:2*ndorb, ndorb:2*ndorb] += soc_d

    emat_n[0:ndorb, 0:ndorb] += soc_d
    emat_n[ndorb:2*ndorb, ndorb:2*ndorb] += soc_d

    # Terms from static core-hole potential
    for i in range(2 * ndorb):
        emat_n[i, i] -= nporb * static_v

    # Crystal field and hoppings between the two Ir-sites
    d = -0.03  # trgional splitting in t2g-subspace 

    # Uncomment the following line to do calculation without hopping and crystal filed splitting.
    t1, t2 = -0.18, 0.036  # hopping between the two-sites in t2g-subspace

    cf_tmp = np.array(
        [ # dzx_1, dzy_1,  dxy_1,    dzx_2, dzy_2,  dxy_2
          [     0,     d,      d,       t1,    t2,     t1],  # dzx_1
          [     d,     0,      d,       t2,    t1,     t1],  # dzy_1
          [     d,     d,      0,       t1,    t1,     t2],  # dxy_1

          [    t1,    t2,     t1,        0,     d,      d],  # dzx_2
          [    t2,    t1,     t1,        d,     0,      d],  # dzy_2
          [    t1,    t1,     t2,        d,     d,      0],  # dxy_2
        ])
    # Including spin degree of freedom, in global axis
    cf_spin = np.zeros((2*ndorb, 2*ndorb), dtype=np.complex)
    cf_spin[0:2*ndorb:2, 0:2*ndorb:2] = cf_tmp
    cf_spin[1:2*ndorb:2, 1:2*ndorb:2] = cf_tmp

    # Transform spin basis to local axis
    # 1/2-spinor matrix
    t_spinor = np.zeros((2*ndorb, 2*ndorb), dtype=np.complex128)
    for i in range(nsite):
        alpha, beta, gamma= edrixs.rmat_to_euler(locaxis[i])
        dmat = edrixs.dmat_spinor(alpha, beta, gamma)
        for j in range(ndorb//2):
            off = i * ndorb + j * 2
            t_spinor[off:off+2, off:off+2] = dmat

    # Transform orbital basis from real cubic to complex harmonics
    t_orb = np.zeros((2*ndorb, 2*ndorb), dtype=np.complex128)
    t_orb[0:ndorb, 0:ndorb] = edrixs.tmat_r2c('t2g', True)
    t_orb[ndorb:2*ndorb, ndorb:2*ndorb] = edrixs.tmat_r2c('t2g', True)
    # Do the tranformation
    cf_spin[:, :] = edrixs.cb_op(cf_spin, np.dot(t_spinor, t_orb))

    emat_i[0:2*ndorb, 0:2*ndorb] += cf_spin
    emat_n[0:2*ndorb, 0:2*ndorb] += cf_spin

    # Write emat and umat to files
    # ED inputs
    edrixs.write_emat(emat_i, "ed/hopping_i.in")
    edrixs.write_umat(umat_i, "ed/coulomb_i.in")

    # XAS inputs
    edrixs.write_emat(emat_n, "xas/hopping_n.in")
    edrixs.write_umat(umat_n, "xas/coulomb_n.in")

    # RIXS inputs
    edrixs.write_emat(emat_i, "rixs_pp/hopping_i.in")
    edrixs.write_umat(umat_i, "rixs_pp/coulomb_i.in")
    edrixs.write_emat(emat_n, "rixs_pp/hopping_n.in")
    edrixs.write_umat(umat_n, "rixs_pp/coulomb_n.in")

    edrixs.write_emat(emat_i, "rixs_ps/hopping_i.in")
    edrixs.write_umat(umat_i, "rixs_ps/coulomb_i.in")
    edrixs.write_emat(emat_n, "rixs_ps/hopping_n.in")
    edrixs.write_umat(umat_n, "rixs_ps/coulomb_n.in")


def get_transop(loc, pos):
    # get dipole transition operator in local axis
    dop = edrixs.get_trans_oper('t2gp32')
    dop_g = np.zeros((2, 3, 6, 4), dtype=np.complex)
    # transform to golobal axis
    for i in range(2):
        for j in range(3):
            for k in range(3):
                dop_g[i, j] += loc[i, j, k] * dop[k]

    # RIXS settings
    thin, thout, phi = 30 / 180.0 * np.pi, 60 / 180.0 * np.pi, 0.0
    ein, eout = 11215.0, 11215.0   # eV
    # Wavevector
    kin, kout = edrixs.get_wavevector_rixs(thin, thout, phi, ein, eout)
    # polarization pi-pi and pi-sigma
    for key, alpha, beta in [('pp', 0, 0), ('ps', 0, np.pi/2.0)]:
        ei, ef = edrixs.dipole_polvec_rixs(thin, thout, phi, alpha, beta)
        T_i = np.zeros((2, 6, 4), dtype=np.complex)
        T_f = np.zeros((2, 4, 6), dtype=np.complex)
        for i in range(2):
            for j in range(3):
                T_i[i] += dop_g[i, j] * ei[j]
                T_f[i] += np.conj(np.transpose(dop_g[i, j] * ef[j]))
            # multiply phase factor
            T_i[i] = T_i[i] * np.exp(+1j * np.dot(pos[i], kin))
            T_f[i] = T_f[i] * np.exp(-1j * np.dot(pos[i], kout))  

        transop_rixs_i = np.zeros((20, 20), dtype=np.complex)
        transop_rixs_f = np.zeros((20, 20), dtype=np.complex)
        for i in range(2):
            off1, off2 = i * 6, 12 + i * 4
            transop_rixs_i[off1:off1 + 6, off2:off2 + 4] = T_i[i]
            transop_rixs_f[off2:off2 + 4, off1:off1 + 6] = T_f[i]
        # write to file
        edrixs.write_emat(transop_rixs_i, "rixs_" + key + "/transop_rixs_i.in")
        edrixs.write_emat(transop_rixs_f, "rixs_" + key + "/transop_rixs_f.in")

    # For XAS, use isotropic polarization
    ei = np.array([1, 1, 1]) / np.sqrt(3.0)
    T_i = np.zeros((2, 6, 4), dtype=np.complex)
    for i in range(2):
        for j in range(3):
            T_i[i] += dop_g[i, j] * ei[j]

    transop_xas = np.zeros((20, 20), dtype=np.complex128)
    for i in range(2):
        off1, off2 = i * 6, 12 + i * 4
        transop_xas[off1:off1+6, off2:off2+4] = T_i[i]
    edrixs.write_emat(transop_xas, "xas/transop_xas.in")

def set_config():
    edrixs.write_config(
        "ed",               # directory to write file
        num_val_orbs=12,    # Number of total valence orbitals,
        num_core_orbs=8,    # Number of total core orbitals
        ed_solver=0,        # Type of ED solver, use arpack
        neval=220,          # Number of eigenvalues obtained,
        nvector=10,         # Number of eigenvectors obtained,
        idump=True,         # Dump eigenvectors to file eigvec.xxx,
        num_gs=10,          # Numer of gound states are used for XAS and RIXS calculations,
        linsys_tol=1E-10,   # Tolerance for the termination of solving linear equations,
        nkryl=500,          # Maximum number of Krylov vectors,
        gamma_in=2.5,       # Core-hole life-time in eV,
        omega_in=-3.9,      # Relative incident x-ray energy at which RIXS calculations are performed,
    )
    shutil.copy('ed/config.in', 'xas/config.in')
    shutil.copy('ed/config.in', 'rixs_pp/config.in')
    shutil.copy('ed/config.in', 'rixs_ps/config.in')


def get_fock_basis(nval_orb, noccu):
    edrixs.write_fock_dec_by_N(nval_orb, noccu, 'ed/fock_i.in')
    shutil.copy('ed/fock_i.in', 'xas/fock_i.in')
    shutil.copy('ed/fock_i.in', 'rixs_pp/fock_i.in')
    shutil.copy('ed/fock_i.in', 'rixs_ps/fock_i.in')
    shutil.copy('ed/fock_i.in', 'rixs_pp/fock_f.in')
    shutil.copy('ed/fock_i.in', 'rixs_ps/fock_f.in')

    edrixs.write_fock_dec_by_N(nval_orb, noccu+1, 'xas/fock_n.in')
    shutil.copy('xas/fock_n.in', 'rixs_pp/fock_n.in')
    shutil.copy('xas/fock_n.in', 'rixs_ps/fock_n.in')


if __name__ == "__main__":
    # site positions, in global xyz-axis
    site_pos = np.array([
                          [2.9164499, 1.6838131, 2.3064262],
                          [2.9164499, 1.6838131, 4.9373740]
                        ])
    # local axis for the two sites, the Wannier orbitals are defined with respect to them
    loc = np.zeros((2, 3, 3), dtype=np.float64)
    loc[0] = np.transpose(np.array([
                                    [+0.70710678, +0.40824829, +0.57735027],   # x
                                    [-0.70710678, +0.40824829, +0.57735027],   # y
                                    [+0.00000000, -0.81649658, +0.57735027]    # z
                                   ]))

    loc[1] = np.transpose(np.array([
                                    [-0.70710678, +0.40824829, -0.57735027],   # x
                                    [+0.70710678, +0.40824829, -0.57735027],   # y
                                    [+0.00000000, -0.81649658, -0.57735027]    # z
                                   ]))

    print("edrixs >>> building directories ...")
    build_dirs()
    print("edrixs >>> Done !")

    print("edrixs >>> set control file config.in ...")
    set_config()
    print("edrixs >>> Done !")

    print("edrixs >>> building fock_basis ...")
    num_val_orbs = 12
    noccu = 9
    get_fock_basis(num_val_orbs, noccu)
    print("edrixs >>> Done !")

    print("edrixs >>> building hopping and coulomb parameters ...")
    get_hopping_coulomb(loc)
    print("edrixs >>> Done !")

    print("edrixs >>> building transition operators from 2p_3/2 to 5d ...")
    get_transop(loc, site_pos)
    print("edrixs >>> Done !")
