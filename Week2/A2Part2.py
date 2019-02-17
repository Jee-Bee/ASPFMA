# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:07:32 2016

@author: Jee-Bee
"""

import numpy as np

def genComplexSine(k, N):
    """
    Inputs:
        k (integer) = frequency index of the complex sinusoid
                      of the DFT
        N (integer) = length of complex sinusoid in samples
    Output:
        The function should return a numpy array
        cSine (numpy array) = The generated complex sinusoid
                              (length N)
    """
    # check input values
    if (N, k >= 0) and (k <= (N-1)):
        # omega = 2 * pi * f == 2 * pi * k
        # t = n * 1 / fs  => -N/2 ... N/2 -1 = 1/ fs
        n = np.arange(0, N)
        cSine = 1 * np.exp(1j * 2 * np.pi * k * n / N)
        return(cSine)
    else:
        print("One of the functions is a non positive value or `k > N-1`")

k = 1
N = 5

complexsinewave = genComplexSine(k, N)
print(complexsinewave)