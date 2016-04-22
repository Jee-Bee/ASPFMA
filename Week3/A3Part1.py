# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 11:42:28 2016

@author: Jee-Bee
"""

import numpy as np
from fractions import gcd
from scipy.fftpack import fft, ifft
import matplotlib.pylab as plt

def minimizeEnergySpreadDFT(x, fs, f1, f2):
    """
    Inputs:
        x (numpy array) = input signal
        fs (float) = sampling frequency in Hz
        f1 (float) = frequency of the first sinusoid component in Hz
        f2 (float) = frequency of the second sinusoid component in Hz
    Output:
        The function should return
        mX (numpy array) = The positive half of the DFT spectrum
                           of the M sample segment of x.
                           mX is (M/2)+1 samples long
                           (M is to be computed)
    """
#    if f1 <= f2:
#        M = fs/(2 * f1)
#    else:
#        M = fs/(2 * f2)
    # M = (f1 * f2)/gcd(f1, f2)
    # @ check calculation description fore M. according to example it have to 
    # be not the right method. but the right method gives wrong value...
    
    M = fs/gcd(f1, f2)
    N = M
    x = x[:M]

    hM1 = np.floor((M+1)/2)
    hM2 = np.floor(M/2)
    dftbuffer = np.zeros(N)
    dftbuffer[:hM1] = x[hM2:]
    dftbuffer[-hM2:] = x[:hM2]

    Mx = fft(dftbuffer)
    Mx = Mx[: len(Mx)/2 + 1]
    return (20 * np.log10(Mx))


# test case 1:
w = 1024.0
fs = 10000.0
f1 = 80.0
f2 = 200.0

t = np.arange(0, w / fs, 1 / fs)
x = np.cos(2 * np.pi * f1 * t) + np.cos(2 * np.pi * f2 * t)
# M = 250 samples
mESX = minimizeEnergySpreadDFT(x, fs, f1, f2)
F = np.arange(0, fs/2, fs / (2 *len(mESX)))
print(len(mESX))

plt.figure()
plt.plot(F, mESX)
plt.show()

# test case 2:
w = 1024.0
fs = 48000.0
f1 = 300.0
f2 = 800.0

t = np.arange(0, w / fs, 1 / fs)
x = np.cos(2 * np.pi * f1 * t) + np.cos(2 * np.pi * f2 * t)
# M = 480 samples
mESX = minimizeEnergySpreadDFT(x, fs, f1, f2)
F = np.arange(0, fs/2, fs / (2 *len(mESX)))
print(len(mESX))

plt.figure()
plt.plot(F, mESX)
plt.show()
