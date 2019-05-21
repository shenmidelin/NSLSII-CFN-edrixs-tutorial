#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

xas = np.loadtxt('xas.dat')

label = [r'$\pi$-pol', r'$\sigma$-pol', r'left-pol', r'right-pol', r'isotropic']

for i in range(len(label)):
    plt.plot(xas[:,0], xas[:,i+1], label=label[i])

plt.xlabel('Incident Energy (eV)')
plt.ylabel('XAS Intensity (a.u.)')

plt.grid()
plt.legend(loc=1, ncol=1, fontsize=12)
plt.subplots_adjust(left=0.15,bottom=0.15)
plt.savefig("xas.pdf")
