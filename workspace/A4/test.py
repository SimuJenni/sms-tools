from A4Part2 import computeSNR
from A4Part3 import computeEngEnv
import os
import sys
import numpy as np
import math
from scipy.signal import get_window
import matplotlib.pyplot as plt

sys.path.append('../../software/models/')
import stft
import utilFunctions as UF
eps = np.finfo(float).eps

inputFile = '../../sounds/piano.wav'
window = 'blackman'
M = 512
N = 1024
H = 128
print computeEngEnv(inputFile, window, 512, 1024, 128)
