#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  5 21:10:50 2018

@author: Jee-Bee
"""

import sys
import os
import numpy as np
# import math
from scipy.signal import get_window
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../sms-tools/software/models/'))
#sys.path.append('~//sms-tools/software/models/')
import stft
import utilFunctions as UF
eps = np.finfo(float).eps


def computeODF(inputFile, window, M, N, H):
    """
    Inputs:
        inputFile (string): input sound file (monophonic with sampling rate of 44100)
        window (string): analysis window type (choice of rectangular, triangular,
            hanning, hamming, blackman, blackmanharris)
        M (integer): analysis window size (odd integer value)
        N (integer): fft size (power of two, such that N > M)
        H (integer): hop size for the STFT computation
    Output:
        The function should return a numpy array with two columns, where the first
        column is the ODF computed on the low frequency band and the second column
        is the ODF computed on the high frequency band.
        ODF[:,0]: ODF computed in band 0 < f < 3000 Hz
        ODF[:,1]: ODF computed in band 3000 < f < 10000 Hz"""
    fs, x = UF.wavread(inputFile)
    w = get_window(window, M)
    if N < M is True:
        raise ValueError("'N' should be greather than 'M'")
    if np.log2(N) % 1 != 0:
        raise ValueError("Input not power of 2")
    Xm, Xp = stft.stftAnal(x, w, N, H)
    Xm = 10 ** (Xm / 20)
    k = Xm.shape[0]
    #f = k × fs / N
    k_1 = np.array([3000, 10000]) * N / fs
    f = fs / 2.0 * np.arange(M) / float(M)
    print(f.shape)
    f_low = np.where(f>3000)[0][0]
    f_high = np.where(f>10000)[0][0]
    print(f_low, f_high, f_high - f_low, k_1)
    engEnv = np.zeros((k, 2))
    engEnv[:,0] = 10 * np.log10(np.sum(np.abs(Xm[:,:f_low]) ** 2, axis=1))
    engEnv[:,1] = 10 * np.log10(np.sum(np.abs(Xm[:,f_low+1:f_high]) ** 2, axis=1))
    engEnv = np.vstack((np.zeros((1,2)), engEnv))
    # ODF = np.zeros((k,2))
    ODF = engEnv[1:-1,:] - engEnv[:-2,:]
    ODF0 = np.where(ODF<0)
    ODF[ODF0[0],ODF0[1]] = 0
    return(ODF)


# Test case 1: 
# Use piano.wav file with
window = 'blackman' 
m = 513
n = 1024  # and 
h = 128  # as input. 
# The bin indexes of the low frequency band span from 1 to 69 (69 samples) and of the high frequency band span from 70 to 232 (163 samples). 
# To numerically compare your output, use loadTestCases.py script to obtain the expected output.
odf = computeODF('../sms-tools/sounds/piano.wav', window, m,n, h)

fs, x = UF.wavread('../sms-tools/sounds/piano.wav')
w = get_window(window, m)
if n < m is True:
    raise ValueError("'N' should be greather than 'M'")
if np.log2(n) % 1 != 0:
    raise ValueError("Input not power of 2")
Xm, Xp = stft.stftAnal(x, w, n, h)
Xm = 10 ** (Xm / 20)

plt.figure(figsize=(10,8))
plt.subplot(2,1,1)
plt.pcolormesh(Xm.T)
plt.xlim(xmin=0, xmax=1350)
#plt.colorbar()
plt.subplot(2,1,2)
plt.plot(odf)
plt.legend(["low freq odf", "high freq odf"])
plt.xlim(xmin=0, xmax=1350)
#plt.ylim(ymin=-10)

#Test case 2:
# Use piano.wav file with 
window = 'blackman'
m = 2047
n = 4096  # and 
h = 128  # as input. 
# The bin indexes of the low frequency band span from 1 to 278 (278 samples) and of the high frequency band span from 279 to 928 (650 samples). 
# To numerically compare your output, use loadTestCases.py script to obtain the expected output.
odf= computeODF('../sms-tools/sounds/piano.wav', window, m,n, h)

fs, x = UF.wavread('../sms-tools/sounds/piano.wav')
w = get_window(window, m)
if n < m is True:
    raise ValueError("'N' should be greather than 'M'")
if np.log2(n) % 1 != 0:
    raise ValueError("Input not power of 2")
Xm, Xp = stft.stftAnal(x, w, n, h)
Xm = 10 ** (Xm / 20)

plt.figure(figsize=(10,8))
plt.subplot(2,1,1)
plt.pcolormesh(Xm.T)
plt.xlim(xmin=0, xmax=1350)
#plt.colorbar()
plt.subplot(2,1,2)
plt.plot(odf)
plt.legend(["low freq odf", "high freq odf"])
plt.xlim(xmin=0, xmax=1350)

# Test case 3:
# Use sax-phrase-short.wav file with 
window = 'hamming'
m = 513
n = 2048  # and 
h = 256  # as input. 
# The bin indexes of the low frequency band span from 1 to 139 (139 samples) and of the high frequency band span from 140 to 464 (325 samples). 
# To numerically compare your output, use loadTestCases.py script to obtain the expected output.
# In addition to comparing results with the expected output, you can also plot your output for these test cases. For test case 1, you can clearly see that the ODFs have sharp peaks at the onset of the piano notes (See (Figure 3). You will notice exactly 6 peaks that are above 10 dB value in the ODF computed on the high frequency band.
odf = computeODF('../sms-tools/sounds/sax-phrase-short.wav', window, m,n, h)

fs, x = UF.wavread('../sms-tools/sounds/sax-phrase-short.wav')
w = get_window(window, m)
if n < m is True:
    raise ValueError("'N' should be greather than 'M'")
if np.log2(n) % 1 != 0:
    raise ValueError("Input not power of 2")
Xm, Xp = stft.stftAnal(x, w, n, h)
Xm = 10 ** (Xm / 20)

plt.figure(figsize=(10,8))
plt.subplot(2,1,1)
plt.pcolormesh(Xm.T)
plt.xlim(xmin=0, xmax=550)
#plt.colorbar()
plt.subplot(2,1,2)
plt.plot(odf)
plt.legend(["low freq odf", "high freq odf"])
plt.xlim(xmin=0, xmax=550)
