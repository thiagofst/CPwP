import numpy as np
import numpy.linalg as la  # Linear Algebra module from numpy
import scipy.linalg as las # Linear Algebra module from SciPy

''' Basic matrix elements '''
# Lets first create a one-dimensional array, that is, a vector
a = np.array([[1,5,-4,-3]])

# To verify the number of dimensions, the shape and also its size:
a.ndim; a.shape; a.size

# Let us now create a two-dimensional array, that is, a matrix
M = np.array([[1,4], [9,2], [-1,6]])

# To select a certain column from M
col = M[:,1]
# Now, to select a certain row
row = M[:1]

# Suppose you have a 1D slice selected from M above. To create a 2D column vector:
column = col.reshape(3,1)

''' Matrix Operations and Functions '''

# Lets define two matrices A and B
A = np.array([[1,2], [1,2]])
B = np.array([[2,1], [2,1]])

# With two matrices, the multiplication can be done as
A @ B

# Matrix Powers. That is, A@A@A@A@A...
la.matrix_power(A, 5)

# Transpose
A.T # Note that A A^t is a symmetric matrix!
# Inverse
A = np.array([[1,2],[3,4]])
# Trace
np.trace(A)
# Determinant
la.det(A)

# Linear systems of equations Ax = b can be solved using: la.solve(A, b)

# Consider the system
# 2x+y−z = 8
# −3x−y+2z− = -11
# 2x+y+2z= -3

# Take the values that multiplies x, y and z and forms a new matrix b with the values on the right hand side of the set of equations
A = np.array([[2, 1, -1], [-3, -1, 2], [-2, 1, 2]])
b = np.array([8, -11, -3])

# The solutions for x, y and z are 2, 3 and -1 respectively

# The Cayley-Hamilton Theorem states that any square matrix satisfies its characteristic polynomial.

A = np.array([[1,2], [3,4]])
trace_A = np.trace(A)
det_A = la.det(A)
I = np.eye(2)
A @ A - trace_A * A + det_A * I

# Pseudo-random matrixes
N = np.random.randint(0,10,[2,2]) # Numeros de 0 a 10 numa matrix 2x2
trace_N = np.trace(N)
det_N = la.det(N)
I = np.eye(2)
N @ N - trace_N * N + det_N * I

def proj(v,w):
	'''Projects vector v into w'''
	v = np.array(v)
	w = np.array(w)
	return (v@w)/(w@w) * w # or np.sum(v * w)/np.sum(w * w) * w
