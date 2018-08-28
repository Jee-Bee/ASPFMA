# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 15:26:57 2016

@author: Jee-Bee
"""

import numpy as np
from scipy.fftpack import fft, ifft
import matplotlib.pyplot as plt


def optimalZeropad(x, fs, f):
    """
    Inputs:
        x (numpy array) = input signal of length M
        fs (float) = sampling frequency in Hz
        f (float) = frequency of the sinusoid in Hz
    Output:
        The function should return
        mX (numpy array) = The positive half of the DFT spectrum
            of the N point DFT after zero-padding input x
            appropriately (zero-padding length to be computed).
            mX is (N/2)+1 samples long
    """
    ps = fs / f  # samples per period
    # M = len(x)
    # Ms = M / ps
    # N = np.ceil(Ms) * ps
    N = np.ceil(len(x) / ps) * ps
    print(N)

    hM1 = np.floor((M+1)/2)
    hM2 = np.floor(M/2)
    dftbuffer = np.zeros(N)
    dftbuffer[:hM1] = x[hM2:]
    dftbuffer[-hM2:] = x[:hM2]

    Mx = fft(dftbuffer)
    Mx = Mx[: len(Mx)/2 + 1]
    return (20 * np.log10(Mx))

# Test case 1:
f = 100.0
M = 25.0
fs = 1000.0
# 5 point zero pad
N = 30

t = np.arange(0, M / fs, 1 / fs)
x = np.cos(2 * np.pi * f * t)
ZEROX = optimalZeropad(x, fs, f)
F = np.arange(0, (fs + 1) / 2, 0.5 * fs / (len(ZEROX) -1))

# bin 3 = 100 Hz
# return is 16 samples
print(len(ZEROX), F[3])

plt.figure()
plt.plot(F, ZEROX)
plt.show()

## Test case 2:
f = 250.0
M = 210.0
fs = 10000.0
# 30 point zero pad
N = 240

t = np.arange(0, M / fs, 1 / fs)
x = np.cos(2 * np.pi * f * t)
ZEROX = optimalZeropad(x, fs, f)
F = np.arange(0, (fs + 1) / 2, 0.5 * fs / (len(ZEROX) -1))

# bin 6 is 250 Hz
# return is 121 samples
print(len(ZEROX), F[6])

plt.figure()
plt.plot(F, ZEROX)
plt.show()
