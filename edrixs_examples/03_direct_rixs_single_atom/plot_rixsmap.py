#!/usr/bin/env python

import numpy as np
import edrixs

rixs_pi = np.loadtxt('rixs_pi.dat')
rixs_sigma = np.loadtxt('rixs_sigma.dat')
eloss = np.loadtxt('eloss.dat')
ominc = np.loadtxt('ominc.dat')

edrixs.plot_rixs_map(rixs_pi, ominc, eloss, "rixs_pi.pdf")
edrixs.plot_rixs_map(rixs_sigma, ominc, eloss, "rixs_sigma.pdf")
