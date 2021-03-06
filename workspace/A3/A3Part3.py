import numpy as np
from scipy.fftpack import fft, fftshift
import math

"""
A3-Part-3: Symmetry properties of the DFT

Write a function to check if the input signal is real and even using the symmetry properties of its DFT. 
The function will return the result of this test, the zerophase windowed version of the input signal 
(dftbuffer), and the DFT of the dftbuffer. 

Given an input signal x of length M, do a zero phase windowing of x without any zero-padding (a dftbuffer, 
on the same lines as the fftbuffer in sms-tools). Then compute the M point DFT of the zero phase windowed 
signal and use the symmetry of the computed DFT to test if the input signal x is real and even. Return 
the result of the test, the dftbuffer computed, and the DFT of the dftbuffer. 

The input argument is a signal x of length M. The output is a tuple with three elements 
(isRealEven, dftbuffer, X), where 'isRealEven' is a boolean variable which is True if x is real 
and even, else False. dftbuffer is the M length zero phase windowed version of x. X is the M point 
DFT of the dftbuffer. 

To make the problem easier, we will use odd length input sequence in this question (M is odd). 

Due to the precision of the FFT computation, the zero values of the DFT are not zero but very small
values < 1e-12 (or -240 dB) in magnitude. For practical purposes, all values with absolute value less 
than 1e-6 (or -120 dB) can be considered to be zero. Use an error tolerance of 1e-6 on the linear 
scale to compare if two floating point arrays are equal. 

EXAMPLE: If x = np.array([ 2, 3, 4, 3, 2 ]), which is a real and even signal, the function returns 
(True, array([ 4., 3., 2., 2., 3.]), array([14.0000+0.j, 2.6180+0.j, 0.3820+0.j, 0.3820+0.j, 2.6180+0.j]))
(values are approximate)

"""
def testRealEven(x):
    """
    Inputs:
        x (numpy array)= input signal of length M (M is odd)
    Output:
        The function should return a tuple (isRealEven, dftbuffer, X)
        isRealEven (boolean) = True if the input x is real and even, and False otherwise
        dftbuffer (numpy array, possibly complex) = The M point zero phase windowed version of x 
        X (numpy array, possibly complex) = The M point DFT of dftbuffer 
    """
    ## Your code here
    
    M = x.size
    h1 = np.floor((M+1)/2)
    h2 = np.floor(M/2)
    x0 = np.zeros(M)
    x0[:h1] = x[h2:]
    x0[h1:] = x[:h2]
    X = fft(x0)
    sym = X[:0:-1]-X[1:]
    isRealEven = (sum(abs(sym))<1e-10) & (sum(abs(np.imag(X)))<1e-10)
    return (isRealEven, x0, X)
