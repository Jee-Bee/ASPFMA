# -*- coding: utf-8 -*-
"""
Created on Wed Jan 02 21:27:37 2018

@author: Jee-Bee
"""

import numpy as np
from scipy.signal import get_window
from scipy.fftpack import fft, fftshift
import math
import matplotlib.pyplot as plt
eps = np.finfo(float).eps

def extractMainLobe(window, M):
    """
    Input:
      window (string): Window type to be used (Either rectangular (‘boxcar’),
          ‘hamming’ or ‘blackmanharris’)
      M (integer): length of the window to be used
    Output:
      The function should return a numpy array containing the main lobe of
      the magnitude spectrum of the window in decibels (dB).
    """
    M = np.asarray(M, int)
    w = get_window(window, M)         # get the window
    N = 8 * M
    W_N = fftshift(fft(w, N))
    W_N = np.abs(W_N)
    W_N[W_N < eps] = eps
    idx_l = int(N / 2)
    lim_l = lim_u = False
    while lim_l is False:
        if W_N[idx_l-1] > W_N[idx_l] and W_N[idx_l+1] > W_N[idx_l]:
            lim_l = True
            print(idx_l)
        else:
            idx_l += -1
    idx_u = int(N / 2)+1
    while lim_u is False:
        if W_N[idx_u-1] > W_N[idx_u] and W_N[idx_u+1] > W_N[idx_u]:
            lim_u = True
            print(idx_u)
        else:
            idx_u += 1
    W_N_dB = 10 * np.log10(W_N)
    # return(W_N_dB, idx_l, idx_u)
    return(W_N_dB[idx_l:idx_u + 1])
    ### Your code here


# Test case 1:
# If you run your code using window = ‘blackmanharris’ and M = 100, the output numpy array should contain 65 samples.
blackmanharris = extractMainLobe('blackmanharris', 100)
# blackmanharris, idx_l, idx_u = extractMainLobe('blackmanharris', 100)
plt.plot(blackmanharris)
# plt.plot(idx_l, blackmanharris[idx_l], 'ro')
# plt.plot(idx_u, blackmanharris[idx_u], 'ro')
print(blackmanharris.shape)
# Test case 2:
# If you run your code using window = ‘boxcar’ and M = 120, the output numpy array should contain 17 samples.
rectangular = extractMainLobe('boxcar', 120)
# rectangular, idx_l, idx_u = extractMainLobe('boxcar', 120)
plt.figure()
plt.plot(rectangular)
# plt.plot(idx_l, rectangular[idx_l], 'ro')
# plt.plot(idx_u, rectangular[idx_u], 'ro')
print(rectangular.shape)
# Test case 3:
# If you run your code using window = ‘hamming’ and M = 256, the output numpy array should contain 33 samples.
hamming = extractMainLobe('hamming', 256)
# hamming, idx_l, idx_u = extractMainLobe('hamming', 256)
plt.figure()
plt.plot(hamming)
# plt.plot(hamming[idx_l:idx_u],'g-')
# plt.plot(idx_l, hamming[idx_l], 'ro')
# plt.plot(idx_u, hamming[idx_u], 'ro')
print(hamming.shape)
