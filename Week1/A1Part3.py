# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 11:24:23 2016

@author: Jee-Bee
"""
#import os
import numpy as np


def hopSamples(x, N):
    """
    Inputs:
        x: input numpy array
        N: a positive integer, (indicating hop size)
    Output:
        A numpy array containing every Nth element in x, starting
        from the first element in x.
    """
    if N > len(x):
        raise ValueError("N is bigger as length of x")
    elif N < 0:
        raise ValueError("N is no positive value")
    elif isinstance(N, np.int) == False:
        raise TypeError("N not of the type int")
    else:
        hop_samples = x[::N]
        return(hop_samples)

#filename = "piano.wav"
#relpath = "sms_tools/sounds/"
#relpath = os.path.normpath(relpath)
#inputfile = os.path.join(os.path.curdir, relpath, filename)

N = 2
x = np.arange(10)

Nth_sample = hopSamples(x, N)
print(Nth_sample)
