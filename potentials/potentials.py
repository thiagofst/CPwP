import matplotlib.pyplot as plt
import numpy as np

plt.rcdefaults(); 
plt.rc('text', usetex = True); 
plt.rc('font', family = 'serif')

# Define two positive and three negative
# charges (x [m], y [m], q [C])
C = [(-3,0,1), (3,0,1), (0,3,-1), (0,-3,-1), (0,-1,-5)] 

[xmin, xmax, ymin, ymax] = [-5, 5, -5, 5]
k = 8.99*10**9  # [Nm^2/C^2], in Coulomb's law

plt.figure(figsize = (8,7))

for i in range(0, len(C)):
    if C[i][2] > 0:
        color = 'bo'
    elif C[i][2] < 0:
        color = 'ro'
    else:
        color = 'wo'
    plt.plot(C[i][0], C[i][1], color)


plt.axis([xmin, xmax, ymin, ymax])
n = 200j  # Mesh grid resolution
Y, X = np.mgrid[xmin:xmax:n, ymin:ymax:n]  # Meshgrid
Ex, Ey = np.array(X*0), np.array(Y*0)

for i in range(0, len(C)):
    R = np.sqrt((X-C[i][0])**2 + (Y-C[i][1])**2)
    Ex = Ex + k*C[i][2]/R**2*(X-C[i][0])/R
    Ey = Ey + k*C[i][2]/R**2*(Y-C[i][1])/R
    

# Plot the result
plt.streamplot(X, Y, Ex, Ey, color='k', density=1, 
               arrowstyle='simple')
plt.xlabel('x [m]', fontsize = 10)
plt.ylabel('y [m]', fontsize = 10)
plt.tight_layout()
plt.savefig('potentials.png', dpi = 200)
plt.show()