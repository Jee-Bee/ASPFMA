#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 07:21:16 2018

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


def computeEngEnv(inputFile, window, M, N, H):
    """
    Inputs:
      inputFile (string): input sound file (monophonic with sampling rate of 44100)
      window (string): analysis window type (choice of rectangular, triangular,
          hanning, hamming, blackman, blackmanharris)
      M (integer): analysis window size (odd positive integer)
      N (integer): FFT size (power of 2, such that N > M)
      H (integer): hop size for the stft computation
    Output:
      The function should return a numpy array engEnv with shape Kx2,
      K = Number of frames
      containing energy envelop of the signal in decibles (dB) scale
      engEnv[:,0]: Energy envelope in band 0 < f < 3000 Hz (in dB)
      engEnv[:,1]: Energy envelope in band 3000 < f < 10000 Hz (in dB)
    """
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
    return(engEnv)



# Example: 
# Running your code on piano.wav file with 
# window = ‘blackman’, 
# M = 513, 
# N = 1024, 
# H = 128, in the plots you can clearly notice the sharp attacks and decay of the piano notes (Figure 2 shows the spectrogram and the energy envelopes).
# In addition, you can also visually analyse which of the two energy envelopes is better for detecting onsets of the piano notes.

#Test case 1:
# Use piano.wav file with
window = 'blackman'
m = 513
n = 1024
h = 128  # as input.
# The bin indexes of the low frequency band span from 1 to 69 (69 samples) and of the high frequency band span from 70 to 232 (163 samples).
# To numerically compare your output, use loadTestCases.py script to obtain the expected output.
engENV= computeEngEnv('../sms-tools/sounds/piano.wav', window, m,n, h)
# print(engENV)

plt.figure(figsize=(8,6))
plt.plot(engENV)


# Test case 2: Use piano.wav file with:
window = 'blackman'
m = 2047
n = 4096
h = 128  # as input. 
# The bin indexes of the low frequency band span from 1 to 278 (278 samples) and of the high frequency band span from 279 to 928 (650 samples).
# To numerically compare your output, use loadTestCases.py script to obtain the expected output.
engENV= computeEngEnv('../sms-tools/sounds/piano.wav', window, m,n, h)
# print(engENV)

plt.figure(figsize=(8,6))
plt.plot(engENV)


# Test case 3:
# Use sax-phrase-short.wav file with
window = 'hamming'
m = 513
n = 2048
h = 256  # as input.
# The bin indexes of the low frequency band span from 1 to 139 (139 samples) and of the high frequency band span from 140 to 464 (325 samples).
# To numerically compare your output, use loadTestCases.py script to obtain the expected output.
engENV= computeEngEnv('../sms-tools/sounds/sax-phrase-short.wav', window, m,n, h)
# print(engENV)

plt.figure(figsize=(8,6))
plt.plot(engENV)
