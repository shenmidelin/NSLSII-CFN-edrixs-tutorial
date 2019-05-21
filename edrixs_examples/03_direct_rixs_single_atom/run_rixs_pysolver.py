#!/usr/bin/env python

import numpy as np
import edrixs

if __name__ == "__main__":
    """
    Ni-3d8: :math:`L_{2,3}`-edge, transition from :math:`2p_{1/2,3/2}\\rightarrow 3d`.

    Use Python ED, XAS, RIXS solvers.

    How to run
    ----------
    python run_rixs_pysolver.py
    """
    # Atomic shell settings
    # ---------------------
    res = edrixs.get_atom_data('Ni', v_name='3d', v_noccu=8, edge='L23')
    name_i, slater_i = [list(i) for i in zip(*res['slater_i'])]  # Initial
    name_n, slater_n = [list(i) for i in zip(*res['slater_n'])]  # Intermediate

    # Slater integrals
    si = edrixs.rescale(slater_i, ([1, 2], [0.65, 0.65]))
    si[0] = edrixs.get_F0('d', si[1], si[2])  # F0_dd

    sn = edrixs.rescale(slater_n, ([1, 2, 4, 5, 6], [0.65, 0.65, 0.95, 0.7, 0.7]))
    sn[0] = edrixs.get_F0('d', sn[1], sn[2])  # F0_dd
    sn[3] = edrixs.get_F0('dp', sn[5], sn[6])  # F0_dp

    slater = (si, sn)

    # Spin-orbit coupling strength
    zeta_d_i, zeta_d_n = res['v_soc_i'][0], res['v_soc_n'][0]
    zeta_p_n = (res['edge_ene'][0] - res['edge_ene'][1]) / 1.5

    # Tetragonal crystal field splitting
    dq10, d1, d3 = 1.3, 0.05, 0.2
    cf = edrixs.cf_tetragonal_d(dq10, d1, d3)

    # Energy shift of the core level
    off = 857.4

    # XAS and RIXS settings
    # ---------------------
    ominc_xas = np.linspace(off - 10, off + 20, 1000)
    thin, thout, phi = 15 / 180.0 * np.pi, 75 / 180.0 * np.pi, 0.0
    poltype_xas = [('linear', 0.0), ('linear', np.pi / 2.0),
                   ('left', 0.0), ('right', 0.0), ('isotropic', 0.0)]
    # L3-edge
    gamma_c, gamma_f = 0.275, 0.1
    ominc_rixs = np.linspace(-5.9 + off, -0.9 + off, 100)
    # L2-edge
    # gamma_c, gamma_f = 0.7, 0.1
    # ominc_rixs = np.linspace(10.9 + off, 14.9 + off, 100)
    eloss = np.linspace(-0.5, 5.0, 1000)
    poltype_rixs = [('linear', 0, 'linear', 0),
                    ('linear', 0, 'linear', np.pi/2.0),
                    ('linear', np.pi/2.0, 'linear', 0.0),
                    ('linear', np.pi/2.0, 'linear', np.pi/2.0)]
    # Run ED
    eval_i, eval_n, trans_op = edrixs.ed_1v1c_py(
        shell_name=('d', 'p'), shell_level=(0, -off), v_soc=(zeta_d_i, zeta_d_n),
        c_soc=zeta_p_n, v_noccu=8, slater=slater, v_cfmat=cf, verbose=1
    )

    # Run XAS
    xas = edrixs.xas_1v1c_py(
        eval_i, eval_n, trans_op, ominc_xas, gamma_c=gamma_c, thin=thin, phi=phi,
        pol_type=poltype_xas, gs_list=[0, 1, 2], temperature=300
    )
    np.savetxt('xas.dat', np.concatenate((np.array([ominc_xas]).T, xas), axis=1))

    # Run RIXS
    rixs = edrixs.rixs_1v1c_py(
        eval_i, eval_n, trans_op, ominc_rixs, eloss, gamma_c=gamma_c, gamma_f=gamma_f, thin=thin,
        thout=thout, phi=phi, pol_type=poltype_rixs, gs_list=[0, 1, 2], temperature=300
    )

    rixs_pi = np.sum(rixs[:, :, 0:2], axis=2)
    rixs_sigma = np.sum(rixs[:, :, 2:4], axis=2)
    np.savetxt('rixs_pi.dat', rixs_pi)
    np.savetxt('rixs_sigma.dat', rixs_sigma)
    np.savetxt('ominc.dat', ominc_rixs)
    np.savetxt('eloss.dat', eloss)
