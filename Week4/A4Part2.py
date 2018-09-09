# -*- coding: utf-8 -*-
"""
Created on Sat Jan 06 09:49:07 2018

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


def computeSNR(inputFile, window, M, N, H):
    """
    Input:
      inputFile (string): input sound file (monophonic with sampling rate of 44100)
      window (string): analysis window type (choice of rectangular, triangular,
          hanning, hamming, blackman, blackmanharris)
      M (integer): analysis window length (odd positive integer)
      N (integer): fft size (power of two, such that N > M)
      H (integer): hop size for the stft computation
    Output:
      The function should return a python tuple of both the SNR values (SNR1, SNR2)
      SNR1 and SNR2 are floats.
    """
    fs, x = UF.wavread(inputFile)
    w = get_window(window, M)
    if N < M is True:
        raise ValueError("'N' should be greather than 'M'")
    if np.log2(N) % 1 != 0:
        raise ValueError("Input not power of 2")
    # Xm, Xp = stft.stftAnal(x, w, N, H)
    # x_ret = stft.stftSynth(Xm, Xp, M, H)
    x_ret = stft.stft(x, w, N, H)
    rest_x_H = x.shape[0] % H
    if rest_x_H % H != 0:
        delta_end = int((H - rest_x_H) // 2)
        print(delta_end)
        x_ret = x_ret[delta_end:-delta_end]
    print(x.shape, x_ret.shape, H, x.shape[0] % H)
    plt.figure(figsize=(8, 6), dpi=120)
    plt.subplot(2,1,1)
    plt.plot(x)
    plt.subplot(2,1,2)
    plt.plot(x_ret)
    E_signal = (np.abs(x) ** 2).sum()
    #E_noise = (np.abs(x - x_ret) ** 2).sum()
    #print(E_signal, E_noise)
    #snr1 = 10 * np.log10(E_signal / E_noise)
    #E_signal = (np.abs(x[M:-M]) ** 2).sum()
    #E_noise = ((np.abs(x - x_ret)[M:-M]) ** 2).sum()
    #print(E_signal, E_noise)
    #snr2 = 10 * np.log10(E_signal / E_noise)
    #return (snr1, snr2)
    return(x, x_ret)

# Test case 1:
# If you run your code using piano.wav file (in the sounds folder of sms-tools) with ‘blackman’ window, M = 513, N = 2048 and H = 128, the output SNR values should be around: (67.57, 304.68)
m = 523
n = 2048
h = 128
path = '.. /sms-tools/sounds/'
file = 'piano.wav'
signal_noise_ratio1, signal_noise_ratio2 = computeSNR('../sms-tools/sounds/piano.wav', 'blackman', m, n, h)
print(signal_noise_ratio1, signal_noise_ratio2)
print('compare with', 67.57, 304.68)
# Test case 2:
# If you run your code using sax-phrase-short.wav file with ‘hamming’ window, M = 512, N = 1024 and H = 64, the output SNR values should be around: (89.51, 306.19)
m = 512
n = 1024
h = 64
file = 'sax-phrase-short.wav'
signal_noise_ratio1, signal_noise_ratio2 = computeSNR('../sms-tools/sounds/sax-phrase-short.wav', 'hamming', m, n, h)
print(signal_noise_ratio1, signal_noise_ratio2)
print('compare with', 89.51, 306.19)
# Test case 3:
# If you run your code using rain.wav file with ‘hann’ window, M = 1024, N = 2048 and H = 128, the output SNR values should be around: (74.63, 304.27)
# Due to precision differences on different machines/hardware, compared to the expected SNR values, your output values can differ by ±10 dB for SNR1 and ±100 dB for SNR2.
m = 1024
n = 2048
h = 128
file = 'rain.wav'
# signal_noise_ratio1, signal_noise_ratio2 = computeSNR('../sms-tools/sounds/rain.wav', 'hann', m, n, h)
print(signal_noise_ratio1, signal_noise_ratio2)
print('compare with', 74.63, 304.27)
