# -*- coding: utf-8 -*-
"""
Created on Tue Apr  5 09:25:13 2016

@author: enjbwink
"""
import os
import sys
# import numpy as np
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), 'sms_tools/software/models/'))
import utilFunctions as uf
#from sms_tools.software.models import utilFunctions
from cython import *

def readAudio(inputFile):
    """
    Input:
        inputFile: the path to the wav file
    Output:
        The function should return a numpy array that
        contains 10 samples of the audio.
    """
#    wavfile = utilFunctions.wavread(inputfile)
    fs, wavfile = uf.wavread(inputfile)
    first_sample = 50001
    wavfile = wavfile[first_sample-1:first_sample + 9]
    return(wavfile)


filename = "piano.wav"
relpath = "sms_tools/sounds/"
relpath = os.path.normpath(relpath)
inputfile = os.path.join(os.path.curdir, relpath, filename)

samples = readAudio(inputfile)
print samples