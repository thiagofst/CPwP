import matplotlib.pyplot as plt
from rungekutta import rungekutta

plt.rcdefaults();
plt.rc('text', usetex = True);
plt.rc('font', family = 'serif')

'''
The Lotka–Volterra equations are a pair of first-order nonlinear differential equations 
used to describe the dynamics of biological systems in which two species interact, one as a predator and the other as prey.'''

# Parameters, RTFM
alpha = 0.1
gamma = 0.05
beta = 0.01
delta = 0.001

# Functions
xdot = lambda t, x, y: alpha*x - x*y*beta
ydot = lambda t, x, y: delta*x*y - gamma*y

# Solver
lotka_volterra = rungekutta(xdot, ydot)
t, y = lotka_volterra.solve([80, 40], .1, 500)

# Results
plt.figure()
plt.grid(True, alpha = 0.4)

plt.plot(t, y[0], 'k', label = r'Prey')
plt.plot(t, y[1], 'r-', label = r'Predator')

plt.xlabel(r'Time/(a.u.)', fontsize = 10)
plt.ylabel(r'Population', fontsize = 10)

plt.legend(fontsize = 10, markerscale = 2, shadow =  'True');
plt.tight_layout(); plt.savefig('lokta_volterra.png', dpi = 200)
plt.close()
