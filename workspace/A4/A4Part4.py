import os
import sys
import numpy as np
from scipy.signal import get_window
import matplotlib.pyplot as plt
import math

sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../software/models/'))
import stft
import utilFunctions as UF
from A4Part3 import computeEngEnv

eps = np.finfo(float).eps

"""
A4-Part-4: Computing onset detection function (Optional)

Write a function to compute a simple onset detection function (ODF) using the STFT. Compute two ODFs one 
for each of the frequency bands, low and high. The low frequency band is the set of all the frequencies 
from 0 - 3000 Hz and the high frequency band is the set of all the frequencies from 3000 - 10000 Hz. 

A brief description of the onset detection function can be found in the pdf document (A4-STFT.pdf, in 
Relevant Concepts section) in the assignment directory (A4). Start with an initial condition of 
ODF(0) = 0 in order to make the length of the ODF same as that of the energy envelope. Remember to 
apply a half wave rectification on the ODF. 

The input arguments to the function are the wav file name including the path (inputFile), window type (window),
window length (M), FFT size (N), and hop size (H). The function should return a numpy array with two columns, 
where the first column is the ODF computed on the low frequency band and the second column is the ODF computed
on the high frequency band.

Use stft.stftAnal() to obtain the STFT magnitude spectrum for all the audio frames. Then compute two 
energy values for each frequency band specified. While calculating frequency bins for each frequency band,
consider only the bins that are within the specified frequency range. For example, for the low frequency 
band consider only the bins with frequency > 0 Hz and < 3000 Hz. This way we also remove the DC offset in 
the signal in energy envelope computation.

To get a better understanding of the energy envelope and its characteristics you can plot the envelopes 
together with the spectrogram of the signal. You can use matplotlib plotting library for this purpose. 
To visualize the spectrogram of a signal, a good option is to use colormesh. You can reuse the code in
sms-tools/lectures/4-STFT/plots-code/spectrogram.py. Either overlay the envelopes on the spectrogram 
or plot them in a different subplot. Make sure you use the same range of the x-axis for both the 
spectrogram and the energy envelopes.

EXAMPLE: Running your code on piano.wav file with window = 'blackman', M = 513, N = 1024, H = 128 in 
the plots you can clearly notice that ODF have sharp peaks at the onset of the piano notes (See figure in 
the accompanying pdf). You will get exactly 6 peaks that are above 10 dB value in the ODF computed on the 
high frequency band. 
"""

def computeODF(inputFile, window, M, N, H):
    """
    Inputs:
            inputFile (string): input sound file (monophonic with sampling rate of 44100)
            window (string): analysis window type (choice of rectangular, triangular, hanning, hamming, 
                blackman, blackmanharris)
            M (integer): analysis window size (odd integer value)
            N (integer): fft size (power of two, bigger or equal than than M)
            H (integer): hop size for the STFT computation
    Output:
            The function should return a numpy array with two columns, where the first column is the ODF 
            computed on the low frequency band and the second column is the ODF computed on the high 
            frequency band.
            ODF[:,0]: ODF computed in band 0 < f < 3000 Hz 
            ODF[:,1]: ODF computed in band 3000 < f < 10000 Hz
    """
    
    ### your code here
    def undoDB(x):
        return np.power(10, np.divide(x,20))
    
    def energy(x, k1, k2):
        x2 = np.power(x[:,k1:k2], 2)
        return np.sum(x2, axis=1)
    
    fs, x = UF.wavread(inputFile)
    w = get_window(window, M)
    bin1 = int(np.ceil(3000*N/fs))
    bin2 = int(np.ceil(10000*N/fs))
    mX, pX = stft.stftAnal(x, fs, w, N, H)
    nrgEnv1 = 10*np.log10( energy(undoDB(mX), 0, bin1))
    nrgEnv2 = 10*np.log10(energy(undoDB(mX), bin1, bin2))    
    engEnv = np.transpose(np.array([nrgEnv1, nrgEnv2]))

    O = engEnv[1:,:]-engEnv[0:-1]
    O[O<0]=0

    return O
    
    """
    plt.subplot(211)
    numFrames = int(mX[:,0].size)
    frmTime = H*np.arange(numFrames-1)/float(fs)                             
    binFreq = np.arange(N/2+1)*float(fs)/N                         
    plt.pcolormesh(frmTime, binFreq, np.transpose(mX[1:,:]))
    plt.title('mX (piano.wav), M=1001, N=1024, H=256')
    plt.autoscale(tight=True)

    plt.subplot(212)
    plt.plot(frmTime, O[:,0], label='low')
    plt.plot(frmTime, O[:,1], label='height')
    plt.title('Energy envelopes')
    plt.autoscale(tight=True)
    
    plt.tight_layout()
    plt.show()
    """
