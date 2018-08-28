# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 10:31:11 2016

@author: Jee-Bee
"""

import os, sys
import numpy as np
from scipy import get_window
# sys.path.append(os.path.dirname(os.path.realpath(__file__)), '../software/models/')
sys.path.append('../sms-tools/software/models')
# import utilFunctions as UF
import dftmodel as DFT


def suppressFreqDFTmodel(x, fs, N):
    """
    Input:
        x (numpy array) = input signal of length M (odd)
        fs (float) = sampling frequency (Hz)
        N (positive integer) = FFT size
    Output:
        The function should return a tuple (y, yfilt)
        y (numpy array) = Output of the dftSynth() without
                        filtering (M samples long)
        yfilt = Output of the dftSynth() with filtering
                (M samples long)
    The first few lines of the code have been written for you,
    do not modify it.
    """
    M = len(x)
    w = get_window('hamming', M)
    outputScaleFactor = w.sum()
    # Make FFT from file
    mX, pX = DFT.dftanal(x, w, N)

    # copy signal before filtering
    mX_filt = mX.copy()
    pX_filt = pX.copy()

    filtLev = 10 ** (- 6)
    filtFreq = 70
    filtFreqBin = np.ceil(filtFreq / (fs / M)) + 1  # calculate frequncy bin higher as 70HZ
    mX_filt[0:filtFreqBin] = filtLev

    # create Inv FFT original signal
    y = DFT.dftSynth(mX, pX, w.size) * outputScaleFactor

    # create Inv FFT filtered signal
    yfilt = DFT.dftSynth(mX_filt, pX_filt, w.size) * outputScaleFactor
    return(y, yfilt)


# Test case 1:
fs = 2048
T = 1
t = np.arange(T * fs) / fs
x = 1 * np.sin(2 * np.pi * 40 * t) + 1 * np.sin(2 * np.pi * 100 * t) + 1 * np.sin(2 * np.pi * 200 * t) + 1 * np.sin(2 * np.pi * 1000 * t)
# yfilt will only contain 100 Hz, 200 Hz and 1000 Hz components.
(y, yfilt) = suppressFreqDFTmodel(x, fs, N)


# Test case 2:
fs = 5000
T = 1
t = np.arange(T * fs) / fs
x = 1 * np.sin(2 * np.pi * 23 * t) + 1 * np.sin(2 * np.pi * 36 * t) + 1 * np.sin(2 * np.pi * 230 * t) + 1 * np.sin(2 * np.pi * 900 * t) + 1 * np.sin(2 * np.pi * 2300 * t)
# yfilt will only contain 230 Hz, 900 Hz and 2300 Hz components.
(y, yfilt) = suppressFreqDFTmodel(x, fs, N)
