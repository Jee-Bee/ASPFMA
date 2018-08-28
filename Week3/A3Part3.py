# -*- coding: utf-8 -*-
"""
Created on Wed Apr 20 15:13:49 2016

@author: Jee-Bee
"""

import numpy as np
from scipy.fftpack import fft
import matplotlib.pyplot as plt


def testRealEven(x):
    """
    Inputs:
        x (numpy array)= input signal of length M (M is odd)
    Output:
        The function should return a tuple (isRealEven, dftbuffer, X)
        isRealEven (boolean) = True if the input x is real and even,
                                and False otherwise
        dftbuffer (numpy array, possibly complex) = The M point zero
                                phase windowed version of x
        X (numpy array, possibly complex) = The M point DFT
                                            of dftbuffer
    """
    N = len(x)
    hM1 = np.floor((N+1)/2)
    hM2 = np.floor(N/2)
    dftbuffer = np.zeros(N)
    dftbuffer[:hM1] = x[hM2:]
    dftbuffer[-hM2:] = x[:hM2]

    Mx = fft(dftbuffer)
    Mx = Mx[: len(Mx)/2 + 1]
    xfh = x[: np.floor(len(x)/2)]  # First half
    # add np.float() otherwise len(x)/ is int so floor
    xsh = x[np.ceil(np.float(len(x))/2):]  # second half
    xsh = xsh[::-1]
    if np.isrealobj(x):
        if all(xfh == xsh):
            return (True, dftbuffer, Mx)
        else:
            return (False, dftbuffer, Mx)
    else:
        return (False, dftbuffer, Mx)

# Test case 1:
x = np.array([2.0, 3.0, 4.0, 3.0, 2.0])
# (True, array([ 4., 3., 2., 2., 3.])
# array([14.0000+0.j, 2.6180+0.j, 0.3820+0.j, 0.3820+0.j, 2.6180+0.j]))
EVENX = testRealEven(x)
print(EVENX)

# Test case 2:
x = np.array([1.0, 2.0, 3.0, 4.0, 1.0, 2.0, 3.0])
# (False, array([ 4., 1., 2., 3., 1., 2., 3.])
# array([ 16.+0.j, 2.+0.69j,2.+3.51j, 2.-1.08j, 2.+1.08j, 2.-3.51j, 2.-0.69j]))
EVENX = testRealEven(x)
print(EVENX)
