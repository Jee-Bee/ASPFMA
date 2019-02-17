# -*- coding: utf-8 -*-
"""
Created on Fri Apr  8 12:07:32 2016

@author: Jee-Bee
"""

import numpy as np

def genSine(A, f, phi, fs, t):
    """
    Inputs:
        A (float) = amplitude of the sinusoid
        f (float) = frequency of the sinusoid in Hz
        phi (float) = initial phase of the sinusoid in radians
        fs (float) = sampling frequency of the sinusoid in Hz
        t (float) = duration of the sinusoid (is second)
    Output:
        The function should return a numpy array
        x (numpy array) = The generated sinusoid (use np.cos())
    """
    # check input values
    if (A, f, fs >= 0) and (f <= (fs/2)):
        t_arr = np.arange(0, t, 1/fs)
        x = A * np.cos(2 * np.pi * f * t_arr + phi)
        return(x)
    else:
        print("One of the functions is a non positive value or `f > fs/2`")

A = 1
f = 10
phi = 1
fs = 50
t = 0.1

sinewave = genSine(A, f, phi, fs, t)
print(sinewave)
