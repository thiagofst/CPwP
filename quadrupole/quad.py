import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D

omega = 1e7 # Frequency
v0 = 5e3 # Velocity of the incident ion
L = 0.1  # [m] Electrodes length
u = 1.66*1e-27  # [m] One atomic unit
q = 1.6022*1e-19  # [C] Elementary charge
r0 = 0.003  # [m] Length from the origin to the electrodes

newparams = {'font.size': 10, 
		  'figure.figsize': (16, 5),
             	  'mathtext.fontset': 'stix', 
             	  'font.family': 'serif',
             	  'lines.linewidth': 2}

plt.rcParams.update(newparams)
plt.rc('text', usetex = True); 
plt.rc('font', family = 'serif')

#
resolution = 200j
equipotential_curves = 30

y, x = np.mgrid[-2*r0:2*r0:resolution, -2*r0:2*r0:resolution]

# Calculating the hyperbolic potential and electric field
Vhyp = (x**2 - y**2)/r0**2
Exhyp = -2/r0**2*x
Eyhyp = 2/r0**2*y

# Calculating the potential and electric field of the point charges
r1 = np.sqrt((x - r0)**2 + y**2)
r2 = np.sqrt((x + r0)**2 + y**2)
r3 = np.sqrt(x**2 + (y - r0)**2)
r4 = np.sqrt(x**2 + (y + r0)**2)
Vlin = -np.log(r1*r2/(r3*r4))
Eylin, Exlin = np.gradient(-Vlin)

plt.figure()
# Plot the hyperbolic potential and electric field
plt.subplot(121)
plt.streamplot(x/r0, y/r0, Exhyp, Eyhyp)
plt.contour(x/r0, y/r0, Vhyp, equipotential_curves)
plt.plot([-1, 1], [0, 0], 'ro')
plt.plot([0, 0], [-1, 1], 'bo')
plt.title('Hyperbolic potential and electric field')
# Plot the potential and electric field of the point charges
plt.subplot(122)
plt.streamplot(x/r0, y/r0, Exlin, Eylin)
plt.contour(x/r0, y/r0, Vlin, equipotential_curves)
plt.plot([-1, 1], [0, 0], 'ro')
plt.plot([0, 0], [-1, 1], 'bo')
plt.title('Potential and electric field of the point charges')
plt.tight_layout()
plt.savefig('fields.png', dpi = 200)
plt.show()