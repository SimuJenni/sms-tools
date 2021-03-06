import numpy as np
from scipy.signal import get_window
from scipy.fftpack import fft
import math
import matplotlib.pyplot as plt
eps = np.finfo(float).eps

""" 
A4-Part-1: Extracting the main lobe of the spectrum of a window

nWrite a function that extracts the main lobe of the magnitude spectrum of a window given a window type 
and its length (M). The function should return the samples corresponding to the main lobe in decibels (dB).

To compute the spectrum, take the FFT size (N) to be 8 times the window length (N = 8*M) 
(For this part, N need not be a power of 2). 

The input arguments to the function are the window type (window) and the length of the window (M). 
The function should return a numpy array containing the samples corresponding to the main lobe of 
the window. 

In the returned numpy array you should include the samples corresponding to both the local minimas
across the main lobe. 

The possible window types that you can expect as input are rectangular ('boxcar'), 'hamming' or
'blackmanharris'.

EXAMPLE: If you run your code using window = 'blackmanharris' and M = 100, the output numpy array 
should contain 65 samples.

NOTE: You can approach this question in two ways: 1) You can write code to find the indices of the 
local minimas across the main lobe. 2) You can manually note down the indices of these local minimas 
by plotting and a visual inspection of the spectrum of the window. If done manually, the indices have 
to be obtained for each possible window types separately (as they differ across different window types).

Tip: log(0) is not well defined, so its a common practice to add a small value such as eps = 1e-16 to the 
magnitude spectrum before computing it in dB. This is optional and will not affect your answers.
"""
def extractMainLobe(window, M):
    """
    Input:
            window (string): Window type to be used (Either rectangular ('boxcar'), 'hamming' or '
                blackmanharris')
            M (integer): length of the window to be used
    Output:
            The function should return a numpy array containing the main lobe of the magnitude spectrum 
            of the window in decibels (dB).
    """

    ### Your code here
    w = get_window(window, M)
    hM1 = int(math.floor((M+1)/2))
    hM2 = int(math.floor(M/2))

    N = 8*M
    fftBuffer = np.zeros(N)
    fftBuffer[:hM1] = w[hM2:]
    fftBuffer[N-hM2:] = w[:hM2]

    X = fft(w, N)
    Xabs = abs(X)
    Xabs[Xabs<eps] = eps
    mX = 20.*np.log10(Xabs)

    hW = min( [i for i in range(len(mX)-1) if mX[i]<mX[i+1] and mX[i]<mX[i-1]])
    result = np.zeros(2*hW + 1)
    result[:hW] = mX[-hW:]
    result[hW:] = mX[:hW+1]
    
    # plt.plot(result)
    # plt.show()

    return result
