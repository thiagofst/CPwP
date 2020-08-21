from matplotlib.pyplot import *
import numpy as np
from scipy.constants import h, c, k
import os

os.chdir('/home/volkoff/VVV/MSC/BlackBody')
rcdefaults(); rc('text', usetex = True); rc('font', family = 'serif')

intensity = lambda wav, T:  (2. * h * c**2) / ( (wav**5) * (np.exp((h * c / (wav * k * T))) - 1.) )
wavs = np.arange(1e-9, 3e-6, 1e-09)

figure()
grid(True, alpha = 0.2); xlim(0, 3000); ylim(0, 2.7e12)
xlabel(r'$\lambda$/nm', fontsize = 10); ylabel(r'B$_{\nu}(\lambda, T)$/W~m$^{-3}$', fontsize = 10)
# Tab. 4.1 - Reid & Hawley - New Light on Dark Stars (2nd Ed.)
intensity2300 = intensity(wavs, 2300.); intensity2400 = intensity(wavs, 2400.); intensity2500 = intensity(wavs, 2500.)
intensity2600 = intensity(wavs, 2600.); intensity2800 = intensity(wavs, 2800.); intensity3100 = intensity(wavs, 3100.)
intensity3259 = intensity(wavs, 3250.); intensity3400 = intensity(wavs, 3400.); intensity3600 = intensity(wavs, 3600.)
intensity3800 = intensity(wavs, 3800.);
#solarintensity = intensity(wavs, 6000.); bsun = solarintensity.max()
# Optical electromagnetic spectrum
axvline(380, lw = 11, alpha = 0.2, color = 'darkmagenta')
axvline(450, lw = 7.5, color = 'navy', alpha = 0.2)
axvline(495, lw = 8, color = 'limegreen', alpha = 0.2)
axvline(570, lw = 12, color = 'yellow', alpha = 0.2)
axvline(590, lw = 8, color = 'darkorange', alpha = 0.2)
axvline(620, lw = 8, color = 'red', alpha = 0.2)
axvline(878, ls = 'dashdot', lw = 1.5, color = 'red')
axvline(2149, ls = 'dashdot', lw = 1.5, color = 'red')
#
plot(wavs*1e+09, intensity2300, 'k-', lw = 1, label = 'Red dwarf stars')
plot(wavs*1e+09, intensity2400, 'k-', lw = 1); plot(wavs*1e+09, intensity2500, 'k-', lw = 1)
plot(wavs*1e+09, intensity2600, 'k-', lw = 1); plot(wavs*1e+09, intensity2800, 'k-', lw = 1)
plot(wavs*1e+09, intensity3100, 'k-', lw = 1); plot(wavs*1e+09, intensity3259, 'k-', lw = 1)
plot(wavs*1e+09, intensity3400, 'k-', lw = 1); plot(wavs*1e+09, intensity3600, 'k-', lw = 1)
#plot(wavs*1e+09, solarintensity, 'b-.', label = r'B$\odot$ = {:e}'.format(bsun))
legend(fontsize = 10, loc = 'upper center', markerscale = 2, shadow = 'True'); tight_layout(); savefig('blackbody.png') #savefig('blackbody_rdw.png', dpi = 300); close()
