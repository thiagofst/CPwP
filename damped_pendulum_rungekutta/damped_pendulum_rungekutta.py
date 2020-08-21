import matplotlib.pyplot as plt
from rungekutta import rungekutta

# Parameters
f  = .1 # Friction
D = 1 # Spring constant
m = 0.1 # Mass of the pendulum

plt.rcdefaults();
plt.rc('text', usetex = True);
plt.rc('font', family = 'serif')

# Equations of motion
ydot = lambda t, x, y: -(f*y + D*x)/m
xdot = lambda t, x, y: -(m*ydot(t, x, y) + D*x)/f

# Solver
oscilator = rungekutta(xdot, ydot)
t, y = oscilator.solve([0, 1], .01, 20)

# Results
plt.figure()
plt.grid(True, alpha = 0.4)

plt.xlabel(r'Time', fontsize = 10)
plt.ylabel(r'Position', fontsize = 10)

plt.plot(t, y[0])
plt.plot(t, y[1])

plt.tight_layout(); plt.savefig('damped_rk.png', dpi = 200)
plt.close()