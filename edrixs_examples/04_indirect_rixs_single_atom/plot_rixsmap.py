#!/usr/bin/env python

import numpy as np
import edrixs

rixs_pi = np.loadtxt('rixs_pi.dat')
eloss = np.loadtxt('eloss.dat')
ominc = np.loadtxt('ominc.dat')

edrixs.plot_rixs_map(rixs_pi, ominc, eloss, "rixs_pi.pdf")
