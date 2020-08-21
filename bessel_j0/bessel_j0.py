import matplotlib.pyplot as plt
import numpy as np
from scipy import special
from mpl_toolkits.mplot3d import Axes3D
import cmocean as oc

plt.rcdefaults(); 
plt.rc('text', usetex = True); 
plt.rc('font', family = 'serif')

def drumhead_height(n, k, distance, angle, t):
   kth_zero = special.jn_zeros(n, k)[-1]
   return np.cos(t) * np.cos(n*angle) * special.jn(n, distance*kth_zero)

theta = np.r_[0:2*np.pi:50j]
radius = np.r_[0:1:50j]
x = np.array([r * np.cos(theta) for r in radius])
y = np.array([r * np.sin(theta) for r in radius])
z = np.array([drumhead_height(1, 1, r, theta, 0.5) for r in radius])

fig = plt.figure(figsize = (8, 7))
ax = Axes3D(fig)
ax.plot_surface(x, y, z, rstride=1, cstride=1, cmap=oc.cm.ice, vmin=-0.5, vmax=0.5)
ax.set_xlabel('X', fontsize = 12)
ax.set_ylabel('Y', fontsize = 12)
ax.set_zlabel('Z', fontsize = 12)
plt.show()