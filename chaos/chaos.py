import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint
from scipy import constants as const

plt.rcdefaults();
plt.rc('text', usetex = True);
plt.rc('font', family = 'serif')

# Define some constants
g = const.g # m/s^2 - Standard acceleration of gravity

def rhs(z, t, l1, l2, m1, m2, g):
    ''' Right hand side of the ODE that describes a Double Pendulum motion '''
    theta1, w1, theta2, w2 = z
    cos12 = np.cos(theta1-theta2)
    sin12 = np.sin(theta1-theta2)
    sin1 = np.sin(theta1)
    sin2 = np.sin(theta2)
    xi = cos12**2 * m2 - m1 - m2
    w1dot = ( L1*m2*cos12*sin12*w1**2 + L2*m2*sin12*w2**2 - m2*g*cos12*sin2  + (m1 + m2)*g*sin1)/(L1*xi)
    w2dot = -( L2*m2*cos12*sin12*w2**2 + L1*(m1 + m2)*sin12*w1**2 + (m1 + m2)*g*sin1*cos12  - (m1 + m2)*g*sin2 )/(L2*xi)
    return w1, w1dot, w2, w2dot

def _cartesian(theta1, w1, theta2, w2, L1, L2):
    ''' Convert theta and omega into cartesian coordinates '''
    x1 = L1 * np.sin(theta1)
    y1 = -L1 * np.cos(theta1)
    x2 = x1 + L2 * np.sin(theta2)
    y2 = y1 - L2 * np.cos(theta2)
    vx1 = L1*np.cos(theta1)*w1
    vy1 = L1*np.sin(theta1)*w1
    vx2 = vx1 + L2*np.cos(theta2)*w2
    vy2 = vy1 + L2*np.sin(theta2)*w2
    return x1, y1, x2, y2, vx1, vy1, vx2, vy2

# Here, we make the (arbitrary) choice 2L1 = L2 and m1 = 3*m2, and we will let the pendulum swing for 50 seconds.

L1, L2 = 1., 2.
m1, m2 = 3., 1.

z0 = [np.pi/2, 0, np.pi/2, 0] # Initial position
tmax, dt = 100, 0.01; t = np.arange(0, tmax+dt, dt)

# Simulation
z = odeint(rhs, z0, t, args = (L1, L2, m1, m2, g))
# Results
theta1, w1, theta2, w2 = z[:,0], z[:,1], z[:,2], z[:,3]
x1, y1, x2, y2, vx1, vy1, vx2, vy2 = _cartesian(theta1, w1, theta2, w2, L1, L2)

# Plot I
plt.figure(); 
plt.grid(True, alpha = 0.4)

plt.xlabel(r'x/L', fontsize = 10)
plt.ylabel(r'y/L', fontsize = 10)

L = 1.1*(L1+L2)
plt.plot(x1, y1, 'k', label = r'm$_{1}$')
plt.plot(x2, y2, 'r', label = r'm$_{2}$')
plt.plot([0, x1[0], x2[0]], [0, y1[0], y2[0]], "-o", label="Initial position")

plt.tight_layout();
plt.legend(fontsize = 10, markerscale = 0.5, shadow =  'True');
plt.savefig('chaos1.png', dpi = 200)

# Plot II
plt.figure(); 
plt.grid(True, alpha = 0.4)

plt.ylabel(r'$\theta$/rad', fontsize = 10)
plt.xlabel(r'$t$/s', fontsize = 10)

plt.plot(t, theta1, 'k', label = r'$\theta_{1}$')
plt.plot(t, theta2, 'r', label = r'$\theta_{2}$')

plt.tight_layout();
plt.legend(fontsize = 10, markerscale = 2, shadow =  'True');
plt.savefig('chaos2.png', dpi = 200)

# Plot III
plt.figure(); 
plt.grid(True, alpha = 0.4)

plt.ylabel(r'$\omega$/(rad/s)', fontsize = 10)
plt.xlabel(r'$t$/s', fontsize = 10)

plt.plot(t, w1, 'k', label = r'$\omega_{1}$')
plt.plot(t, w2, 'r', label = r'$\omega_{2}$')

plt.tight_layout();
plt.legend(fontsize = 10, markerscale = 2, shadow =  'True');
plt.savefig('chaos3.png', dpi = 200)

# Plot IIII
# Phase-space diagram for the double pendulum
plt.figure(); 
plt.grid(True, alpha = 0.4)

plt.xlabel(r'$\theta_{i}$/rad', fontsize = 10)
plt.ylabel(r'$\omega_{i}$/(rad/s)', fontsize = 10)

plt.title(r"Phase-space diagram, $\theta_{10}=%.1f$, $\theta_{20}=%.1f$ "%(theta1[0], theta2[0]) + r"$\omega_{10}=%.1f$, $\omega_{20}=%.1f$"%(w1[0], w2[0]), loc = 'left')
plt.plot(theta1, w1, 'k-', label = r'$i=1$')
plt.plot(theta2, w2, 'r-', label = r'$i=2$')

plt.tight_layout();
plt.legend(fontsize = 10, markerscale = 2, shadow =  'True');
plt.savefig('chaos4.png', dpi = 200)

plt.close()