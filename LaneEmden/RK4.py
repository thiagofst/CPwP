import numpy as np
import matplotlib.pyplot as plt
import os

plt.rcdefaults(); plt.rc('text', usetex = True); plt.rc('font', family = 'serif')

'''
It solves a differential equation like y' - f(x,y) = 0 at some point xf
For the Lane-Emden equation, we should first makethe substitution Z = d\xi/d\theta
Then we get
z' = -Dn^n - 2z/x
'''

laneemden = lambda xi,Dn,z,n: -1.*np.power(Dn, n)-(2./xi)*z # Lane-Emden equation (LE)

def rk4(func, h, x0, y0, z0, xf):
	''' 4th order Runge-Kutta method (RK4)'''
	yn = y0
	xn = x0
	zn = z0
	yvals = []
	xvals = []
	while xn < xf:
		k1 = h*func(xn, yn, zn)
		k2 = h*func(xn+0.5*h, yn+0.5*k1, zn)
		k3 = h*func(xn+0.5*h, yn+0.5*k2, zn)
		k4 = h*func(xn+h, yn+k3, zn)
		zn += (k1 + k2 + k3 + k4)/6.0
		yn += h*zn
		xn += h
		yvals.append(yn)
		xvals.append(xn)
	return (xvals, yvals)

'''
It's well known that LE is analytical just for n = 0, 1 and 5
Let's make a comparisson between analytical and numerical solutions
'''

xi = np.linspace(0, 10, 100) # Choose a set of adimensional radius

# Analytical solutions
# Chandrasekhar (1967), An introduction to the study of stellar structure

le0 = lambda xi: 1. - ((1./6.)*xi**2.) 			# n = 0
le1 = lambda xi: np.sin(xi)/xi 					# n = 1
le5 = lambda xi: 1./np.sqrt(1. + ((xi**2)/3)) 	# n = 5

theta0 = le0(xi) # theta_0 (xi)
theta1 = le1(xi) # theta_1 (xi)
theta5 = le5(xi) # theta_5 (xi)

# Numerical solutions for the ''analytical indexes''' via RK4 method
lane0 = lambda x,y,z: laneemden(x, y, z, 0); (x_, y_) = rk4(lane0, 1e-5, 1-((1e-10)/6), 1, 0, 10)
lane1 = lambda x,y,z: laneemden(x, y, z, 1); (x1, y1) = rk4(lane1, 1e-5, 1-((1e-10)/6), 1, 0, 10)
lane5 = lambda x,y,z: laneemden(x, y, z, 5); (x5, y5) = rk4(lane5, 1e-5, 1-((1e-10)/6), 1, 0, 10)

# Plot the results
linestyle_tuple = [('densely dashdotted', (0, (3, 1, 1, 1))),
				   ('dashdotted', (0, (3, 5, 1, 5))),] # Just some line_styles

plt.figure(); plt.grid(True, alpha = 0.5); plt.xlim(0, 10); plt.ylim(-0.5, 1.)
plt.xlabel(r'$\xi$', fontsize = 10)
plt.ylabel(r'$\theta_{n}(\xi)$', fontsize = 10)
# n = 0
plt.plot(x_, y_, ls = 'dashdot', lw = 1., color = 'black', label = r'$n$ = 0')
plt.plot(xi, theta0, ls = 'dashdot', lw = 1., color = 'black', label = '__nolabel__')
# n = 1
plt.plot(x1, y1, ls = (0, (3,1,1,1)), lw = 1., color = 'black', label = r'$n$ = 1')
plt.plot(xi, theta1, ls = (0, (3,1,1,1)), lw = 1., color = 'black', label = '__nolabel__')
#n = 5
plt.plot(x5, y5, ls = (0, (3,5,1,5)), lw = 1., color = 'black', label = r'$n$ = 5')
plt.plot(xi, theta5, ls = (0, (3,5,1,5)), lw = 1., color = 'black', label = '__nolabel__')
#
plt.legend(fontsize = 10, loc = 'upper right', markerscale = 2, shadow = 'True')
plt.tight_layout(); plt.savefig('analytical_vs_numerical.png', dpi = 200); plt.show()

# Just the analytical solutions
plt.figure(); plt.grid(True, alpha = 0.5); plt.xlim(0, 10); plt.ylim(-0.5, 1.)
plt.xlabel(r'$\xi$', fontsize = 10)
plt.ylabel(r'$\theta_{n}(\xi)$', fontsize = 10)
# n = 0
plt.plot(xi, theta0, ls = 'dashdot', lw = 1., color = 'black', label = r'$n$ = 0')
# n = 1
plt.plot(xi, theta1, ls = (0, (3,1,1,1)), lw = 1., color = 'black', label = r'$n$ = 1')
#n = 5
plt.plot(xi, theta5, ls = (0, (3,5,1,5)), lw = 1., color = 'black', label = r'$n$ = 5')
#
plt.legend(fontsize = 10, loc = 'upper right', markerscale = 2, shadow = 'True')
plt.tight_layout(); plt.savefig('analytical.png', dpi = 200); plt.close()


# Numerical solutions for non-integrers polytropes of order n = 0.5, 1.5, 2.5,..., 4.5
lane05 = lambda x,y,z: laneemden(x, y, z, 0.5); (x_5, y_5) = rk4(lane05, 1e-5, 1-((1e-10)/6), 1, 0, 6)
lane15 = lambda x,y,z: laneemden(x, y, z, 1.5); (x15, y15) = rk4(lane15, 1e-5, 1-((1e-10)/6), 1, 0, 6) # Compressed stars (such as neutron stars) in their final stages of evolution
lane25 = lambda x,y,z: laneemden(x, y, z, 2.5); (x25, y25) = rk4(lane25, 1e-5, 1-((1e-10)/6), 1, 0, 6)
lane35 = lambda x,y,z: laneemden(x, y, z, 3.5); (x35, y35) = rk4(lane35, 1e-5, 1-((1e-10)/6), 1, 0, 6)
lane45 = lambda x,y,z: laneemden(x, y, z, 4.5); (x45, y45) = rk4(lane45, 1e-5, 1-((1e-10)/6), 1, 0, 6)
lane55 = lambda x,y,z: laneemden(x, y, z, 5.5); (x55, y55) = rk4(lane55, 1e-5, 1-((1e-10)/6), 1, 0, 6)
plt.figure(); plt.grid(True, alpha = 0.4)
plt.xlim(1, 5); plt.ylim(0, 1.)
plt.xlabel(r'$\xi$', fontsize = 10)
plt.ylabel(r'$\theta_{n}(\xi)$', fontsize = 10)
#
plt.plot(x_5, y_5, lw = 1.5, c = 'black', label = r'$n$ = 0.5')
plt.plot(x15, y15, lw = 1.5, c = 'red', label = r'$n$ = 1.5')
plt.plot(x25, y25, lw = 1.5, c = 'yellow', label = r'$n$ = 2.5')
plt.plot(x35, y35, lw = 1.5, c = 'limegreen', label = r'$n$ = 3.5')
plt.plot(x45, y45, lw = 1.5, c = 'navy', label = r'$n$ = 4.5')
plt.plot(x55, y55, lw = 1.5, c = 'darkmagenta', label = r'$n$ = 4.5')
#
plt.legend(fontsize = 10, loc = 'upper right', markerscale = 2, shadow = 'True')
plt.tight_layout(); plt.savefig('rk4_non-integers.png', dpi = 100); plt.close()


# Numerical solutions
lane0 = lambda x,y,z: laneemden(x, y, z, 0.); (x_, y_) = rk4(lane0, 1e-5, 1-((1e-10)/6), 1, 0, 10)
lane1 = lambda x,y,z: laneemden(x, y, z, 1.); (x1, y1) = rk4(lane1, 1e-5, 1-((1e-10)/6), 1, 0, 10) # Compressed stars (such as neutron stars) in their final stages of evolution
lane2 = lambda x,y,z: laneemden(x, y, z, 2.); (x2, y2) = rk4(lane2, 1e-5, 1-((1e-10)/6), 1, 0, 10)
lane3 = lambda x,y,z: laneemden(x, y, z, 3.); (x3, y3) = rk4(lane3, 1e-5, 1-((1e-10)/6), 1, 0, 10)
lane4 = lambda x,y,z: laneemden(x, y, z, 4.); (x4, y4) = rk4(lane4, 1e-5, 1-((1e-10)/6), 1, 0, 10)
lane5 = lambda x,y,z: laneemden(x, y, z, 5.); (x5, y5) = rk4(lane5, 1e-5, 1-((1e-10)/6), 1, 0, 10)

plt.figure(figsize = (12, 4)); plt.grid(True, alpha = 0.4)
plt.xlim(1, 10); plt.ylim(-1.5, 1.)
plt.xlabel(r'$\xi$', fontsize = 12)
plt.ylabel(r'$\theta_{n}(\xi)$', fontsize = 12)
#
plt.plot(x_, y_, lw = 1., c = 'black', ls = 'dashdot', label = r'$n$ = 0')
plt.plot(x1, y1, lw = 1., c = 'red', ls = 'dashdot', label = r'$n$ = 1')
plt.plot(x2, y2, lw = 1., c = 'yellow', ls = 'dashdot', label = r'$n$ = 2')
plt.plot(x3, y3, lw = 1., c = 'limegreen', ls = 'dashdot', label = r'$n$ = 3')
plt.plot(x4, y4, lw = 1., c = 'navy', ls = 'dashdot', label = r'$n$ = 4')
plt.plot(x5, y5, lw = 1., c = 'darkmagenta', ls = 'dashdot', label = r'$n$ = 4')
#
plt.legend(fontsize = 10, loc = 'lower left', markerscale = 2, shadow = 'True')
plt.tight_layout(); plt.savefig('rk4.png', dpi = 200); plt.close()
