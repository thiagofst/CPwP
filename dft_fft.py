import numpy as np
from numpy import random as rnd
from scipy import fftpack as fft

#import timeit # Just to measure the execution time of the code
import matplotlib.pyplot as plt

def DFT(x):
    ''' One-dimensional Discrete Fourier Transform (DFT) of a vector 'x'
    Parameters:
    x: double array
    Returns:
    y: double array. The Fourier transform of 'x'
    '''
    n = len(x); y = [0]*n
    omega = np.exp(-2.0j*np.pi/n)
    for k in range(0, n):
        y[k] = np.sum(x*omega**(np.arange(0, n)*k))
    return y


def inverseDFT(y):
    ''' One-dimensional inverse DFT of 'y'
    Parameters:
    y: double array
    Returns:
    x: The inverse FT of 'y'
    '''
    n = len(y); x = [0]*n
    omega = np.exp(2.j*np.pi/n)
    for k in range(0, n):
        x[k] = np.sum(y*omega**(np.arange(0, n)*k))/float(n)
    return x

# Cooley-Tuckey algorithm

def CooleyTukeyRadix2FFT(x):
    ''' One-dimensional Fast Fourier Transform (FFT) of 'x' using radix-2 Cooley-Tuckey (CT) algorithm.
        In this case, the vector to be transformed must have a power of 2 number of elements
    Parameters:
    x: double array
    Returns a recursive formula for calculating the FFT
    '''
    if ( len(x) & (len(x) - 1)):
        raise Exception("The number of elements in x has to be a power of 2!")
    # Recursive formula for calculating the FFT.
    def foo(x):
        n = len(x)
        if n == 1:
            y = x
        else:
            y2 = foo(x[0:n:2])
            y1 = foo(x[1:n + 1:2])
            d = np.exp(-2j*np.pi/n)**np.arange(0,n/2)
            y = np.append(y2 + d*y1,y2 - d*y1)
        return y
    return foo(x)

def inverseCooleyTukeyRadix2FFT(y):
    ''' One-dimensional inverse FFT of 'x' using radix-2 CT algorithm. '''
    if (len(y) & (len(y) - 1)):
        raise Exception("The number of elements in x has to be a power of 2!")
    def foo(y):
        n = len(y)
        if n == 1:
            x = y
        else:
            x2 = foo(y[0:n:2])
            x1 = foo(y[1:n + 1:2])
            d = np.exp(2j*np.pi/n)**np.arange(0,n/2)
            x = np.append(x2 + d*x1,x2 - d*x1)
        return x
    return foo(y)/len(y)


x = rnd.randint(500, size = 2**15)
y = CooleyTukeyRadix2FFT(x)
x_ = inverseCooleyTukeyRadix2FFT(y)
