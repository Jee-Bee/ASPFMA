# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 11:16:41 2016

@author: Jee-Bee
"""
import os
import sys
import numpy as np
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'sms_tools/software/models/'))
import utilFunctions as uf


def minMaxAudio(inputFile):
    """
    Input:
        inputFile: file path to the wav file
    Output:
        A tuple of the minimum and the maximum value of the audio
        samples, like: (min_val, max_val)
    """
    fs, wavfile = uf.wavread(inputFile)
    minval = np.min(wavfile)
    maxval = np.max(wavfile)
    return((minval, maxval))


filename = "oboe-A4.wav"
relpath = "sms_tools/sounds/"
relpath = os.path.normpath(relpath)
inputfile = os.path.join(os.path.curdir, relpath, filename)

minmax = minMaxAudio(inputfile)
print minmax
