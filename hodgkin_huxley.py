import matplotlib.pyplot as plt
from math import *
from rungekutta import rungekutta

plt.rcdefaults();
plt.rc('text', usetex = True);
plt.rc('font', family = 'serif')


'''
The Hodgkin–Huxley model is a set of nonlinear differential equations that approximates the electrical characteristics of excitable cells such as neurons and cardiac myocytes.
Basically, It describes how action potentials in neurons are initiated and propagated. 
Alan Hodgkin and Andrew Huxley received the 1963 Nobel Prize in Physiology or Medicine for this work.

The lipid bilayer is represented as a capacitance (Cm)
Voltage-gated ion channels are represented by electrical conductances (gn, where n is the specific ion channel) that depend on both voltage and time. 
Leak channels are represented by linear conductances (gL). 
The electrochemical gradients driving the flow of ions are represented by voltage sources (En) whose voltages are determined by the ratio of the intra- and extracellular concentrations of the ionic species of interest. 
Finally, ion pumps are represented by current sources (Ip).
The membrane potential is denoted by Vm.

I is the total membrane current per unit area, Cm is the membrane capacitance per unit area, gK and gNa are the potassium and sodium conductances per unit area, respectively, VK and VNa are the potassium and sodium reversal potentials, respectively, and gl and Vl are the leak conductance per unit area and leak reversal potential, respectively. The time dependent elements of this equation are Vm, gNa, and gK, where the last two conductances depend explicitly on voltage as well.

'''

# Parameters, RTFM
C = 1;
V = 0;
n = 0.32;
m = 0.08;
h = 0.6;
I = 30;
gK = 36;
gNa = 120;
gL = 0.3;
EK = -12;
ENa = 120;
EL = 10.6;
c = 0

def I(t):
	global c
	if ((round(t/10)+10)%3)  >=2 and c<=10:
		c =+ 0.01
		return 15
	else:
		if not ((round(t/10)+10)%3)  >=2:
			c = 0
		return  0

# Functions
alpha_n  = lambda V: .01 * (10-V)/(exp(1 - .1 * V) -1)
alpha_m = lambda V: .1 * (25 - V)/(exp(2.5 - 0.1*V) - 1)
alpha_h = lambda V: .07 * exp(-V/20)

beta_n = lambda V: 0.125 * exp(-V / 80)
beta_m = lambda V: 4 * exp(-V / 18)
beta_h = lambda V: 1 / (exp(3 - 0.1 * V) + 1)

vdot = lambda t, v, n, m, h: (I(t) - gK*(n**4)*(v - EK) - gNa*(m**3)*h*(v - ENa) - gL*(v - EL))/C
ndot = lambda t, v, n, m, h: alpha_n(v) * (1 - n) - beta_n(V) * n
mdot = lambda t, v, n, m, h: alpha_m(v) * (1 - m) - beta_m(V) * m
hdot = lambda t, v, n, m, h: alpha_h(v) * (1 - h) - beta_h(V) * h

hh = rungekutta(vdot, ndot, mdot, hdot)
t, y = hh.solve([V, n, m, h], .01, 50)

# Results
plt.figure(figsize = (8,7))
plt.subplot(311)
plt.grid(True, alpha = 0.4)
plt.plot(t, y[0])
plt.ylabel(r'Membrane potential/mV', fontsize = 10)

plt.subplot(312)
plt.grid(True, alpha = 0.4)
plt.plot(t, y[0])
plt.plot(t, y[2])
plt.plot(t, y[3])
plt.ylabel(r'Gating', fontsize = 10)

plt.subplot(313)
plt.grid(True, alpha = 0.4)
plt.plot(t, [I(k) for k in t])
plt.ylabel(r'External current', fontsize = 10)
plt.xlabel(r'Time', fontsize = 10)
plt.tight_layout()
plt.savefig('hodgkin_huxley.png', dpi = 200)
