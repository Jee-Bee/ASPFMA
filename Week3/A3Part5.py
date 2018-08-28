# -*- coding: utf-8 -*-
"""
Created on Thu Apr 21 10:41:51 2016

@author: Jee-Bee
"""

import os, sys
import numpy as np
import matplotlib.pyplot as plt
from scipy import get_window
sys.path.append(os.path.dirname(os.path.realpath(__file__)), '../software/models/')
# import utilFunctions as UF
import dftmodel as DFT

def zpFFTsizeExpt(x, fs):
    """
    Inputs:
        x (numpy array) = input signal (2*M = 512 samples long)
        fs (float) = sampling frequency in Hz
    Output:
        The function should return a tuple (mX1_80, mX2_80, mX3_80)
        mX1_80 (numpy array): The first 80 samples of the magnitude
                                spectrum output of dftAnal for Case-1
        mX2_80 (numpy array): The first 80 samples of the magnitude
                                spectrum output of dftAnal for Case-2
        mX3_80 (numpy array): The first 80 samples of the magnitude
                                spectrum output of dftAnal for Case-3

    The first few lines of the code to generate xseg and the windows
    have been written for you, please use it and do not modify it.
    """
    M = len(x)/2
    xseg = x[:M]
    w1 = get_window('hamming',M)
    w2 = get_window('hamming',2*M)

    mX1_80 = Mx1[:80]
    mX2_80 = Mx2[:80]
    mX3_80 = Mx3[:80]
    return(mX1_80, mX2_80, mX3_80)

# Test case 1:
# x (of length 512 samples),
# output is a tuple where mX1 80, mX2 80, mX3 80 are the 
rst 80 samples of the magnitude spectrum
N = 512 # length samples
fs = 1000 # Hz
f = 110 # Hz
t = np.arange(0,N,)

