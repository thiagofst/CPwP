import matplotlib.pyplot as plt
import numpy as np
import astropy.constants
import scipy.constants as const
import os
from astropy.io import ascii
from scipy.constants import Rydberg as rh

plt.rcdefaults();
plt.rc('text', usetex = True);
plt.rc('font', family = 'serif') # LaTeX style

os.chdir('/home/volkoff/VVV/MSC/Astro')

plt.figure(figsize = (12, 7))
plt.subplot(211)

s = ascii.read('/home/volkoff/VVV/MSC/Astro/sun-like.dat')
(wav, flux) = (s['lk'], s['flux'])

balmer = lambda n: 1./(rh*(1/4 - 1/n**2))
n = np.arange(3,10,1)
lines = balmer(n)/(1e-10)

plt.grid(True, alpha = 0.4);
plt.title(r'G2$-$V star', fontsize = 12, loc = 'left')
plt.xlim(3500, 8000);
plt.ylim(0.2, 1.2)
plt.xlabel(r'$\lambda$/$\AA$', fontsize = 12)
plt.ylabel(r'Flux/(erg/(cm$^{2}$ s $\AA$)', fontsize = 12)
plt.plot(wav, flux, 'k-', label = '__nolabel__')
#plt.axvline(3800, lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'black', label = '__nolabel__')
#plt.axvline(7500, lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'black', label = '__nolabel__')
plt.axvline(lines[0], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'red', label = r'H$_{\alpha}$')
plt.axvline(lines[1], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'orange', label = r'H$_{\beta}$')
plt.axvline(lines[2], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'yellow', label = r'H$_{\gamma}$')
plt.axvline(lines[3], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'limegreen', label = r'H$_{\delta}$')
plt.axvline(lines[4], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'cyan', label = r'H$_{\varepsilon}$')
plt.axvline(lines[5], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'navy', label = r'H$_{\zeta}$')
plt.axvline(lines[6], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'darkmagenta', label = r'H$_{\eta}$')
plt.legend(fontsize = 10, markerscale = 2, loc = 'lower right', shadow = 'True')

plt.subplot(212)

s = ascii.read('/home/volkoff/VVV/MSC/Astro/vega-like.dat')
(wav, flux) = (s['lk'], s['flux'])

balmer = lambda n: 1./(rh*(1/4 - 1/n**2))
n = np.arange(3,10,1)
lines = balmer(n)/(1e-10)

plt.grid(True, alpha = 0.4);
plt.title(r'A0$-$V star', fontsize = 12, loc = 'left')
plt.xlim(3500, 8000);
#plt.ylim(0, 2.6)
plt.xlabel(r'$\lambda$/$\AA$', fontsize = 12)
plt.ylabel(r'Flux/(erg/(cm$^{2}$ s $\AA$)', fontsize = 12)
plt.plot(wav, flux, 'k-', label = '__nolabel__')
#plt.axvline(3800, lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'black', label = '__nolabel__')
#plt.axvline(7500, lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'black', label = '__nolabel__')
plt.axvline(lines[0], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'red', label = r'H$_{\alpha}$')
plt.axvline(lines[1], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'orange', label = r'H$_{\beta}$')
plt.axvline(lines[2], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'yellow', label = r'H$_{\gamma}$')
plt.axvline(lines[3], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'limegreen', label = r'H$_{\delta}$')
plt.axvline(lines[4], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'cyan', label = r'H$_{\varepsilon}$')
plt.axvline(lines[5], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'navy', label = r'H$_{\zeta}$')
plt.axvline(lines[6], lw = 1.5, alpha = 0.5, ls = 'dashdot', color = 'darkmagenta', label = r'H$_{\eta}$')
plt.legend(fontsize = 10, markerscale = 2, loc = 'lower right', shadow = 'True')

plt.tight_layout(); plt.savefig('vega-sun-like.eps'); plt.show()
