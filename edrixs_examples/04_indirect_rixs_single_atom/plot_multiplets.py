#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

eigvals = np.loadtxt('eigvals.dat')

plt.ylim((-0.1,3))
plt.grid()
plt.plot(eigvals[:,0], eigvals[:,1] - min(eigvals[:,1]), '-o', linewidth=1, markersize=3)
plt.xlabel(r'Multiplets')
plt.ylabel(r'Energy (eV)')
plt.savefig('multiplets.pdf')
