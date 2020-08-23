import matplotlib.pyplot as plt
import numpy as np
import os
import numpy.fft
from scipy import signal

os.chdir('/home/volkoff/fourier')

plt.rcdefaults(); 
plt.rc('text', usetex = True); 
plt.rc('font', family = 'serif')

# Set
x = np.arange(-10, 10, 1e-5)

# Gaussian distribution
gaussian = lambda x: np.exp(-1*(x)**2/2)
fft_gaussian = np.fft.fftshift(np.abs(np.fft.fft(gaussian(x)))) / np.sqrt(len(gaussian(x)))

plt.figure(figsize = (12, 4))
plt.subplot(121) # Beta numbers normally distributed
plt.title('Normal distribution', fontsize = 12, loc = 'center')
plt.grid(True, alpha = 0.4)
plt.xlabel(r'$\beta$', fontsize = 12)
plt.ylabel(r'$\mathcal{N}(\beta)$', fontsize = 12)
plt.plot(x,  gaussian(x), 'k-', lw = 1.5)
plt.subplot(122) # FFT
plt.grid(True, alpha = 0.4)
plt.xlim(-0.002, 0.002)
plt.xlabel(r'$\beta$', fontsize = 12)
plt.ylabel(r'$\mathcal{F}\{\mathcal{N}(\beta)\}$', fontsize = 12)
fft_gaussian = np.fft.fftshift(np.abs(np.fft.fft(gaussian(x)))) / np.sqrt(len(gaussian(x)))
plt.plot(x, fft_gaussian, 'k-', lw = 1.5)
plt.tight_layout()
plt.savefig('gaussian.png', dpi  = 200)
plt.close() # plt.show()

# Rectangular function
rect = lambda x: np.where(abs(x)<=0.5, 1, 0)
y = np.zeros(len(x)); 
y[200:400] = 1
ffty = np.fft.fft(y); ffty = np.fft.fftshift(ffty)
fft_rect = np.real(ffty)

plt.figure(figsize = (12, 4))
plt.subplot(121) # Beta numbers 
plt.grid(True, alpha = 0.4)
plt.xlim(-2.5, 2.5)
plt.title('Rectangular function', fontsize = 12, loc = 'center')
plt.xlabel(r'$\beta$', fontsize = 12)
plt.ylabel(r'$\Pi$($\beta$)', fontsize = 12)
plt.plot(x,  rect(x), 'k-', lw = 1.5)
plt.subplot(122) # FFT
plt.grid(True, alpha = 0.4)
plt.xlim(-.3, .3)
plt.xlabel(r'$\beta$', fontsize = 12)
plt.ylabel(r'$\mathcal{F}\{\Pi(\beta)\}$', fontsize = 12)
plt.plot(x, fft_rect, 'k-', lw = 1.5)
plt.tight_layout()
plt.savefig('rectangular.png', dpi  = 200)
plt.close() # plt.show()

# Triangular function
triang = lambda x: signal.triang(x)
x_ = triang(100)
A = numpy.fft.fft(x_, 2048) / (len(x_)/2.0)
freq = np.linspace(-0.5, 0.5, len(A))
response = 20 * np.log10(np.abs(np.fft.fftshift(A / np.abs(A).max())))

plt.figure(figsize = (12, 4))
plt.subplot(121) # Beta numbers 
plt.grid(True, alpha = 0.4)
plt.title('Triangular function', fontsize = 12, loc = 'center')
plt.xlabel(r'$\beta$', fontsize = 12)
plt.ylabel(r'$\Lambda$($\beta$)', fontsize = 12)
plt.plot(x_, 'k-', lw = 1.5)
plt.subplot(122) # FFT
plt.grid(True, alpha = 0.4)
plt.ylim(-140, 10)
plt.xlim(-.3, .3)
plt.xlabel(r'$\beta$', fontsize = 12)
plt.ylabel(r'$\mathcal{F}\{\Lambda(\beta)\}$', fontsize = 12)
plt.plot(freq, response, 'k-', lw = 1.5)
plt.tight_layout()
plt.savefig('triangular.png', dpi  = 200)
plt.close() # plt.show()

# Sawtooth function
saw = lambda x: signal.sawtooth(x)
saw_fft_ = np.fft.fft(saw(x)); saw_fft = np.real(saw_fft_)

plt.figure(figsize = (12, 4))
plt.subplot(121) # Beta numbers 
plt.grid(True, alpha = 0.4)
plt.title('Sawtooth function', fontsize = 12, loc = 'center')
plt.xlabel(r'$\beta$', fontsize = 12)
plt.ylabel(r'$\mathcal{S}$($\beta$)', fontsize = 12)
plt.plot(x, saw(x), 'k-', lw = 1.5)
plt.subplot(122) # FFT
plt.grid(True, alpha = 0.4)
plt.xlabel(r'$\beta$', fontsize = 12)
plt.ylabel(r'$\mathcal{F}\{\mathcal{S}(\beta)\}$', fontsize = 12)
plt.plot(x, saw_fft, response, 'k-', lw = 1.5)
plt.tight_layout()
plt.savefig('sawtooth.png', dpi  = 200)
plt.close() # plt.show()

# Shah function
plt.rcdefaults(); 
sha = [u'\u0428']

Nx, dx = 200, 10 # Shah function
x1 = []
y1 = []

for i in range(-Nx, Nx+1):
    x1.append(i)
    if np.mod(np.abs(i), dx)<1:
        y1.append(1.)
    else:
        y1.append(0.)

xshah, yshah = np.array(x1), np.array(y1)

nftpoints = 100
F = numpy.fft.fft(y1, norm = 'ortho', n=nftpoints)
freq = numpy.fft.fftshift(numpy.fft.fftfreq(nftpoints))

plt.figure(figsize = (12, 4))
plt.subplot(121); 
plt.grid(True, alpha = 0.4);
plt.title('Dirac comb function', fontsize = 12, loc = 'center')
plt.xlabel('X', fontsize = 12)
plt.ylabel('%s(x)'%(sha[0]), fontsize = 12)
plt.plot(x1, y1, 'k-', lw = 1.)
plt.subplot(122); 
plt.grid(True, alpha = 0.4);
plt.xlabel('X', fontsize = 12)
plt.ylabel('FT{%s(x)}'%(sha[0]), fontsize = 12)
plt.plot(freq, np.absolute(F), 'k-', lw = 1.)
plt.tight_layout();
plt.savefig('dirac.png', dpi  = 200)
plt.close() # plt.show()

#

plt.figure(figsize = (6, 7))

plt.rcdefaults(); 

sha = [u'\u0428']

Nx, dx = 200, 10 # Shah function
x1 = []
y1 = []

for i in range(-Nx, Nx+1):
    x1.append(i)
    if np.mod(np.abs(i), dx)<1:
        y1.append(1.)
    else:
        y1.append(0.)

xshah, yshah = np.array(x1), np.array(y1)

nftpoints = 100
F = numpy.fft.fft(y1, norm = 'ortho', n=nftpoints)
freq = numpy.fft.fftshift(numpy.fft.fftfreq(nftpoints))

plt.subplot(421)
plt.grid(True, alpha = 0.4);
plt.title('Dirac comb function', fontsize = 10, loc = 'center')
plt.xlabel('X', fontsize = 10)
plt.ylabel('%s(x)'%(sha[0]), fontsize = 10)
plt.plot(x1, y1, 'k-', lw = 1.)

plt.subplot(422)
plt.grid(True, alpha = 0.4);
plt.xlabel('X', fontsize = 10)
plt.ylabel('FT{%s(x)}'%(sha[0]), fontsize = 10)
plt.plot(freq, np.absolute(F), 'k-', lw = 1.)

plt.rcdefaults(); 
plt.rc('text', usetex = True); 
plt.rc('font', family = 'serif')

plt.subplot(423)

plt.title('Normal distribution', fontsize = 10, loc = 'center')
plt.grid(True, alpha = 0.4)
plt.xlabel(r'$\beta$', fontsize = 10)
plt.ylabel(r'$\mathcal{N}(\beta)$', fontsize = 10)
plt.plot(x,  gaussian(x), 'k-', lw = 1.5)

plt.subplot(424)
plt.grid(True, alpha = 0.4)
plt.xlim(-0.002, 0.002)
plt.xlabel(r'$\beta$', fontsize = 10)
plt.ylabel(r'$\mathcal{F}\{\mathcal{N}(\beta)\}$', fontsize = 10)
fft_gaussian = np.fft.fftshift(np.abs(np.fft.fft(gaussian(x)))) / np.sqrt(len(gaussian(x)))
plt.plot(x, fft_gaussian, 'k-', lw = 1.5)

plt.subplot(425)
plt.grid(True, alpha = 0.4)
plt.xlim(-2.5, 2.5)
plt.title('Rectangular function', fontsize = 10, loc = 'center')
plt.xlabel(r'$\beta$', fontsize = 10)
plt.ylabel(r'$\Pi$($\beta$)', fontsize = 10)
plt.plot(x,  rect(x), 'k-', lw = 1.5)

plt.subplot(426)
plt.grid(True, alpha = 0.4)
plt.xlim(-.3, .3)
plt.xlabel(r'$\beta$', fontsize = 10)
plt.ylabel(r'$\mathcal{F}\{\Pi(\beta)\}$', fontsize = 10)
plt.plot(x, fft_rect, 'k-', lw = 1.5)

plt.subplot(427)
plt.grid(True, alpha = 0.4)
plt.title('Triangular function', fontsize = 10, loc = 'center')
plt.xlabel(r'$\beta$', fontsize = 10)
plt.ylabel(r'$\Lambda$($\beta$)', fontsize = 10)
plt.plot(x_, 'k-', lw = 1.5)

triang = lambda x: signal.triang(x)
x_ = triang(100)
A = numpy.fft.fft(x_, 2048) / (len(x_)/2.0)
freq = np.linspace(-0.5, 0.5, len(A))
response = 20 * np.log10(np.abs(np.fft.fftshift(A / np.abs(A).max())))

plt.subplot(428)
plt.grid(True, alpha = 0.4)
plt.ylim(-140, 10)
plt.xlim(-.3, .3)
plt.xlabel(r'$\beta$', fontsize = 10)
plt.ylabel(r'$\mathcal{F}\{\Lambda(\beta)\}$', fontsize = 10)
plt.plot(freq, response, 'k-', lw = 1.5)

plt.tight_layout()
plt.savefig('fourier.png', dpi  = 200)
plt.show()