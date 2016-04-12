# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:07:32 2016

@author: Jee-Bee
"""

import numpy as np

def IDFT(X):
    """
    Input:
        X (numpy array) = frequency spectrum (length N)
    Output:
        The function should return a numpy array of length N
        x (numpy array) = The IDFT of the frequency spectrum X
                          (length N)
    """
    N = len(X)
    k = np.arange(N)
    print(k)
    x = []
    for n in range(N):
        s = np.exp(1j * 2 * np.pi * n / N * k)
        x = np.append(x, 1 / N * np.sum(X * s))
    return(x)

X = np.array([1, 1, 1, 1])
x = IDFT(X)
print(x)
