#!/usr/bin/env python

import numpy as np
import collections
import edrixs
from mpi4py import MPI

if __name__ == "__main__":
    '''
    U-5f2 :math:`L_{3}`-edge, :math:`2p_{3/2}\\rightarrow 6d` transition, indirect RIXS.

    Orbital order: 5f, 6d, 2p32
    '''
    # mpi4py env
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    size = comm.Get_size()

    # Occupancy of U 5f orbitals
    noccu = 2
    res = edrixs.get_atom_data('U', v_name=('5f', '6d'), v_noccu=(noccu, 0),
                               edge='L3', trans_to_which=2, label=('f', 'd', 'p'))
    if rank == 0:
        print(res, flush=True)

    # Slater integrals
    si = collections.OrderedDict(res['slater_i'])
    sn = collections.OrderedDict(res['slater_n'])

    # Initial Hamiltonian, 5f-5f
    si['F2_ff'] = si['F2_ff'] * 0.77
    si['F4_ff'] = si['F4_ff'] * 0.77
    si['F6_ff'] = si['F6_ff'] * 0.77
    si['F0_ff'] = edrixs.get_F0('f', si['F2_ff'], si['F4_ff'], si['F6_ff'])

    # Intermediate Hamiltonian, 5f-5f
    sn['F2_ff'] = sn['F2_ff'] * 0.77
    sn['F4_ff'] = sn['F4_ff'] * 0.77
    sn['F6_ff'] = sn['F6_ff'] * 0.77
    sn['F0_ff'] = edrixs.get_F0('f',  sn['F2_ff'], sn['F4_ff'], sn['F6_ff'])
    # 5f-6d
    sn['F0_fd'] = edrixs.get_F0('fd', sn['G1_fd'], sn['G3_fd'], sn['G5_fd'])
    # 5f-2p
    sn['F0_fp'] = edrixs.get_F0('fp', sn['G2_fp'], sn['G4_fp'])
    # 6d-2p
    sn['F0_dp'] = edrixs.get_F0('dp', sn['G1_dp'], sn['G3_dp'])

    slater = (list(si.values()), list(sn.values()))

    # Spin-Orbit Coupling (SOC) zeta
    # 5f
    zeta_f_i = res['v_soc_i'][0] * 0.9
    zeta_f_n = res['v_soc_n'][0] * 0.9
    # 6d
    zeta_d_i = res['v_soc_i'][1]
    zeta_d_n = res['v_soc_n'][1]

    # RIXS settings
    thin, thout, phi = 45 / 180.0 * np.pi, 45 / 180.0 * np.pi, 0.0
    gamma_c = res['gamma_c'][0]
    gamma_f = 0.1
    ominc_xas = np.linspace(-10, 20, 1000)
    ominc_rixs = np.linspace(5, 15, 10)
    eloss = np.linspace(-0.2, 3.0, 1000)

    poltype_xas = [('isotropic', 0.0)]

    # pi-pi and pi-sigma
    poltype_rixs = [('linear', 0.0, 'linear', 0.0),
                    ('linear', 0.0, 'linear', np.pi / 2.0)]

    shell_name = ('f', 'd', 'p32')

    # Run ED
    eval_i, denmat = edrixs.ed_2v1c_fort(
        comm, shell_name, shell_level=(0, 10.0, 0),
        v1_soc=(zeta_f_i, zeta_f_n), v2_soc=(zeta_d_i, zeta_d_n),
        v_tot_noccu=noccu, slater=slater, ed_solver=2, neval=100,
        nvector=9, ncv=300, idump=True
    )

    # Run XAS
    xas, xas_poles = edrixs.xas_2v1c_fort(
        comm, shell_name, ominc_xas, gamma_c=gamma_c,
        v_tot_noccu=noccu, trans_to_which=2, thin=thin, phi=phi,
        pol_type=poltype_xas, num_gs=9, nkryl=500, temperature=300
    )
    if rank == 0:
        np.savetxt('xas.dat', np.concatenate((np.array([ominc_xas]).T, xas), axis=1))
        edrixs.dump_poles(xas_poles, 'xas_poles')

    # Run RIXS
    rixs, rixs_poles = edrixs.rixs_2v1c_fort(
        comm, shell_name, ominc_rixs, eloss, gamma_c=gamma_c, gamma_f=gamma_f,
        v_tot_noccu=noccu, trans_to_which=2, thin=thin, thout=thout, phi=phi,
        pol_type=poltype_rixs, num_gs=9, nkryl=500, linsys_max=1000, temperature=300
    )

    if rank == 0:
        edrixs.dump_poles(rixs_poles, 'rixs_poles')
        rixs_pi = np.sum(rixs[:, :, 0:2], axis=2)
        np.savetxt('rixs_pi.dat', rixs_pi)
        np.savetxt('ominc.dat', ominc_rixs)
        np.savetxt('eloss.dat', eloss)
