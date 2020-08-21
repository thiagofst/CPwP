import matplotlib.pyplot as plt
import numpy as np
import astropy.constants
import scipy.constants as const
import os
from astropy.io import ascii

plt.rcdefaults();
plt.rc('text', usetex = True);
plt.rc('font', family = 'serif') # LaTeX style

os.chdir('/home/volkoff/VVV/MSC/Astro')

saha = lambda zi,zii,chi,t,ne: (1/ne) * ((alpha * t)**1.5) * ((2 * zii) / zi) * (np.exp(-chi/(k*t))) # Saha equation

# Initial guess for electron density
ne = (lambda atoms, vol: (atoms/vol)/500)(10e20,1.3e4)
# Constant
alpha = 1.8e10 #2 * pi * electron mass * k / h^2
# Partition functions
zii = 1.; zi = 2.
# The energy of the first excited state of hydrogen
chi = 13.6 # eV

k =  .00008617385 #const.physical_constants['Boltzmann constant in eV/K'][0]
t = np.arange(3e3, 25e3) # Temperatures

r = saha(zi,zii,chi,t,ne)
rr = r/(1+r)
rrr = r/(1+r) * 1/(1+r)

# Nii/Ntotal for Hydrogen from the Saha equation
plt.figure(); plt.grid(True, alpha = 0.4); plt.xlim(5000, 22.5e3)
plt.ylabel(r'[N$_{II}$/N$_{total}]_{H}$', fontsize = 12);
plt.xlabel(r'Temperature/K', fontsize = 12)
plt.plot(t, rr, 'k-', lw = 1.3, ls = 'solid')
#
plt.axvline(8.3e3, lw = 1.3, ls = 'dotted', color = 'black', label = r'5$\%$ of ionised  H$_{2}$')
plt.axvline(9.6e3, lw = 1.3, ls = 'dashdot', color = 'black', label = r'50$\%$ of ionised  H$_{2}$')
plt.axvline(11.3e3, lw = 1.3, ls = 'dashed', color = 'black', label = r'95$\%$ of ionised  H$_{2}$')
plt.axvline(9520, lw = 1.3, ls = 'solid', color = 'red', label = 'Maximum intensity of Balmer lines')
#
plt.legend(fontsize = 10, markerscale = 2, loc = 'lower right', shadow = 'True')
plt.tight_layout(); plt.savefig('saha.eps'); plt.show()

#
plt.figure(); plt.grid(True, alpha = 0.4); plt.xlim(6000, 18000)
plt.ylabel(r'[N$_{2}$/N$_{total}]_{H}$', fontsize = 12);
plt.xlabel(r'Temperature/K', fontsize = 12)
plt.plot(t, rrr, 'k-', lw = 1.3, ls = 'solid')
#
plt.axvline(t[np.argmax(rrr)], lw = 1.3, ls = 'dashdot', color = 'red', label = 'Hydrogen atom ionisation temperature \n T = %s K'%(t[np.argmax(rrr)]))
#
plt.legend(fontsize = 10, markerscale = 2, loc = 'upper right', shadow = 'True')
plt.tight_layout(); plt.savefig('saha-II.eps'); plt.show()

