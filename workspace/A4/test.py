from A4Part2 import computeSNR

inputFile = '../../sounds/piano.wav'
window = 'blackman'
M = 512
N = 1024
H = 128
print computeSNR(inputFile, window, 512, 1024, 128)
