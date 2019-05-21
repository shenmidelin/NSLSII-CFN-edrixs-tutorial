#!/usr/bin/env python

import numpy as np
import matplotlib.pyplot as plt

eval_i = np.loadtxt('eval_i.dat')

plt.grid()
plt.plot(eval_i[:,0], eval_i[:,1] - min(eval_i[:,1]), '-o')
plt.xlabel(r'Multiplets')
plt.ylabel(r'Energy (eV)')
plt.savefig('multiplets.pdf')
