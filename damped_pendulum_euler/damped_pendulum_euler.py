# Damped pendulum: 'A real pendulum cannot oscilate forever fue to friction!'

import matplotlib.pyplot as plt
import numpy as np
from scipy import constants as const

plt.rcdefaults();
plt.rc('text', usetex = True);
plt.rc('font', family = 'serif')

# Define some constants
g = const.g # m/s^2 - Standard acceleration of gravity

def rhs(theta, w, dt):
    ''' Gives the RHS of the ODE that describes a damped simple pendulum '''
    dw = -np.sin(theta)*dt*g/L - b/(L*m)*w*dt
    dtheta = w*dt
    return dtheta, dw

def euler_step(theta,w,dt):
    ''' One step of the Euler method '''
    dtheta, dw = rhs(theta, w, dt)
    w = w + dw
    theta = theta + dtheta
    return theta, w

def euler(theta0, w0, dt, n):
    ''' Solution via Euler method '''
    theta = (n+1)*[0]
    w = (n+1)*[0]
    #
    theta[0] = theta0
    w[0] = w0
    for k in range(n):
        theta[k+1], w[k+1] = euler_step(theta[k], w[k], dt)
    return theta, w

# Define some constants
m = 1.      # kg. Mass of rod
L = 1.      # m. Length of rod
w0 = 10     # 1/s. Initial angular velocity
theta0 = 3. # rad. Initial launching angle
T = 20.     # s. Time of simulation
n = 100000  # Step number
b = .5      # kg m. Damping factor

t = np.linspace(0, T, n + 1)
theta, _ = euler(theta0, w0, 0.0002, n)

#
plt.figure();
plt.grid(True, alpha = 0.4)

plt.title(r'Angular position of the rod', fontsize = 12)
plt.xlabel(r'$t$/s', fontsize = 12)
plt.ylabel(r'$\theta$/rad', fontsize = 12)

plt.plot(t, theta, 'k', lw = 1.5, ls = 'solid')
plt.tight_layout(); plt.savefig('angular_rod.png', dpi = 200)

# Conservation of Energy

U = lambda theta: m * g * L * (1 - np.cos(theta)) # Potential energy
K = lambda w: 0.5 * m * L ** 2*np.array(w)**2 # Kinectic energy

#
plt.figure()
plt.grid(True, alpha = 0.4)

plt.xlim(0, 12.5)
plt.ylim(0, 40)

plt.xlabel(r'$t$/s', fontsize = 10)
plt.ylabel(r'$E$/J', fontsize = 10)

plt.plot(t, U(theta), label = r'Potential energy')
plt.plot(t, K(_), label = 'Kinectic energy')
plt.plot(t, (U(theta)+K(_)), label = 'Total energy')

plt.legend(fontsize = 10, markerscale = 2, shadow =  'True', loc = 'center');
plt.tight_layout(); plt.savefig('energy.png', dpi = 200)
