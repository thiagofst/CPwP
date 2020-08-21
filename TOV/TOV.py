# https://colab.research.google.com/drive/1yMD2j3Y6TcsykCI59YWiW9WAMW-SPf12#scrollTo=MQsQ6z8RgqOp

import matplotlib.pyplot as plt
import numpy as np

plt.rcdefaults(); plt.rc('text', usetex = True); plt.rc('font', family = 'serif')

lcgs = 1.476701332464468e+05

def EOS_p2erc(p, K = 100, Gamma = 2.):
    ''' Given pressure, returns energy density, rest-mass density and sound speed '''
    ene = (p/K)**(1./Gamma) + p/(Gamma-1.)
    rho = (p/K)**(1./Gamma)
    cs2 = K*Gamma*(Gamma-1)/(Gamma-1 + K*Gamma*rho**(Gamma-1))*rho**(Gamma-1)
    return ene, rho, cs2


def EOS_r2pe(rho, K = 100., Gamma = 2.):
  ''' Given rest-mass density, return energy density and pressure
   Polytropic EOS: P = k rho^Gamma '''
  p = K*rho**Gamma
  e = rho + p/(Gamma-1.);
  return p, e


def TOV(t, y):
  ''' Tolmann-Oppenheimer-Volkhoff (TOV) equations
  d/dt y(t) = R.H.S.
  '''
  r = t
  #
  m = y[0] # mass of a sphere of radius r
  p = y[1] # pressure
  #
  ene, dummy1, dummy2 = EOS_p2erc(p)
  #
  dy = np.empty_like(y)
  dy[0] = 4*np.pi*ene*r**2
  dy[1] = -(ene+p)*(m + 4*np.pi*r**3*p)/(r*(r-2*m))
  return dy


def found_radius(t,y ,pfloor=0.):
  '''
  Event function: Zero of pressure
  ODE integration stops when this function returns True
  '''
  return ((y[1]-pfloor)<=0.)

def solve_ode_euler(t, y0, dydt_fun, stop_event=None, verbose=False):
  """
  Euler algorithm
  NOTE: solution is not stored/saved, just return last point
  """
  N = len(t)
  dt = np.diff(t)[0] # Assume a uniformly space-t array
  y = y0
  for i in range(N):
    yprev = np.copy(y) # Store previous for returning pre-event data
    y += dt * dydt_fun(t[i],y)
    if verbose: print(t[i],y)
    if stop_event:
      if bool(stop_event(t[i],y)):
        print('Event reached')
        return t[i-1], yprev
  if stop_event: print('No event reached')
  return t[i], y

rmin, rmax = 1e-6, 20.
rspan = np.linspace(rmin, rmax, 100)

rho0    = 1.28e-3 # Central (maximal) rest-mass density
p0,e0   = EOS_r2pe(rho0)
m0      = 4./3.*np.pi*e0*rmin**3
sol0    = [m0, p0]

t, sol = solve_ode_euler(rspan, sol0, TOV, stop_event=found_radius, verbose=True)

# Mass and radius
R = t * lcgs * 1e-5 # Km
M = sol[0] # Msun
pmin = sol[1]
print(pmin,R,M)

def set_initial_conditions(rho, rmin):
  ''' Utility routine to set initial data, given rho_0
  https://arxiv.org/pdf/1305.3510.pdf
  https://stellarcollapse.org/nsmasses
  '''
  p,e = EOS_r2pe(rho)
  m    = 4./3.*np.pi*e*rmin**3
  return m, p


rspan = np.linspace(rmin, rmax, 1000)
rhospan = np.linspace(0.6e-4,7e-3,200)
R = []
M = []

for rho0 in rhospan:
    sol0 = set_initial_conditions(rho0, rmin)
    t, sol = solve_ode_euler(rspan, sol0, TOV, stop_event=found_radius)
    R.append(t)
    M.append(sol[0])

M = np.array(M)
R = np.array(R)

immax = np.argmax(M)
Mmax = M[immax]
Rmax = R[immax]
print(immax,Mmax,Rmax)

km = lcgs * 1e-5

plt.figure(); plt.grid(True, alpha = 0.4)
plt.title(r'Mass$-$Radius diagram for Neutron Stars', fontsize = 10, loc = 'right')
plt.xlabel('$R$ $(km)$', fontsize = 10)
plt.ylabel('$M$ $(M_\odot)$', fontsize = 10)
plt.plot(R*km,M, 'k--', alpha = 0.3, label = '__nolabel__')
plt.plot(Rmax*km,Mmax, 'ro', label = r'Maximum NS mass')
plt.legend(fontsize = 10, loc = 'upper right', markerscale = 0.8, shadow = 'True')
plt.tight_layout(); plt.savefig('tov.png', dpi = 200); plt.close()
