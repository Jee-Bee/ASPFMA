# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:07:32 2016

@author: Jee-Bee
"""

import numpy as np

def DFT(x):
    """
    Input:
        x (numpy array) = input sequence of length N
    Output:
        The function should return a numpy array of length N
        X (numpy array) = The N point DFT of the input sequence x
    """
    N = len(x)
    n = np.arange(N)
    X = []
    for k in range(N):
        s = np.exp(1j * 2 * np.pi * k * n / N)
        X = np.append(X, np.sum(x * np.conjugate(s)))
    return(X)

x = np.array([1, 2, 3, 4])
X = DFT(x)
print(X)
