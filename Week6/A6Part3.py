import os
import sys
import numpy as np
import math
from scipy.signal import get_window
import matplotlib.pyplot as plt

# sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../../software/models/'))
sys.path.append(os.path.join(os.path.dirname(os.path.realpath(__file__)), '../sms-tools/software/models/'))
import utilFunctions as UF
import harmonicModel as HM
import stft

eps = np.finfo(float).eps

"""
A6Part3 - Compute amount of inharmonicity present in a sound

Write a function that measures the amount of inharmonicity present in a pitched/harmonic sound. The 
function should measure the mean inharmonicity in the sound over the time interval t1 to t2.

The input argument to the function are the wav file name including the path (inputFile), start (t1) 
and end time (t2) of the audio segment to compute inharmonicity, analysis window (window), window 
size (M), FFT size (N), hop size (H), error threshold used in the f0 detection (f0et), magnitude 
threshold for spectral peak picking (t), minimum allowed f0 (minf0), maximum allowed f0 (maxf0) and 
number of harmonics to be considered in the computation of inharmonicity (nH). The function returns 
a single numpy float, which is the mean inharmonicity over time t1 to t2. 

A brief description of the method to compute inharmonicity is provided in the Relevant Concepts 
section of the assignment pdf. The steps to be done are:
1. Use harmonicModelAnal function in harmonicModel module for computing the harmonic frequencies and 
their magnitudes at each audio frame. The first harmonic is the fundamental frequency. For 
harmonicModelAnal use harmDevSlope=0.01, minSineDur=0.0. Use harmonicModelAnal to estimate harmonic 
frequencies and magnitudes for the entire audio signal.
2. For the computation of the inharmonicity choose the frames that are between the time interval 
t1 and t2. Do not slice the audio signal between the time interval t1 and t2 before estimating 
harmonic frequencies. 
3. Use the formula given in the Relevant Concepts section to compute the inharmonicity measure for the 
given interval. Note that for some frames some of the harmonics might not be detected due to their low 
energy. For handling such cases use only the detected harmonics (and set the value of R in the equation 
to the number of detected hamonics) to compute the inharmonicity measure. 
All the detected harmonics have a non-zero frequency.

In this question we will work with a piano sound ('../../sounds/piano.wav'), a typical example of an 
instrument that exhibits inharmonicity 
(http://en.wikipedia.org/wiki/Piano_acoustics#Inharmonicity_and_piano_size). 

Test case 1: If you run your code with inputFile = '../../sounds/piano.wav', t1=0.2, t2=0.4, window='hamming', M=2047, N=2048, H=128, f0et=5.0, t=-90, minf0=130, maxf0=180, nH = 25, the returned output should be 1.4543. 

Test case 2: If you run your code with inputFile = '../../sounds/piano.wav', t1=2.3, t2=2.55, window='hamming', M=2047, N=2048, H=128, f0et=5.0, t=-90, minf0=230, maxf0=290, nH = 15, the returned output should be 1.4874. 

Test case 3: If you run your code with inputFile = '../../sounds/piano.wav', t1=2.55, t2=2.8, window='hamming', M=2047, N=2048, H=128, f0et=5.0, t=-90, minf0=230, maxf0=290, nH = 5, the returned output should be 0.1748. 

Optional/Additional tasks
An interesting task would be to compare the inharmonicities present in the sounds of different instruments. 
"""
#def estimateInharmonicity(inputFile = '../../sounds/piano.wav', t1=0.1, t2=0.5, window='hamming',
def estimateInharmonicity(inputFile = '../sms-tools/sounds/piano.wav', t1=0.1, t2=0.5, window='hamming',
                            M=2048, N=2048, H=128, f0et=5.0, t=-90, minf0=130, maxf0=180, nH = 10):
    """
    Function to estimate the extent of inharmonicity present in a sound
    Input:
        inputFile (string): wav file including the path
        t1 (float): start time of the segment considered for computing inharmonicity
        t2 (float): end time of the segment considered for computing inharmonicity
        window (string): analysis window
        M (integer): window size used for computing f0 contour
        N (integer): FFT size used for computing f0 contour
        H (integer): Hop size used for computing f0 contour
        f0et (float): error threshold used for the f0 computation
        t (float): magnitude threshold in dB used in spectral peak picking
        minf0 (float): minimum fundamental frequency in Hz
        maxf0 (float): maximum fundamental frequency in Hz
        nH (integer): number of integers considered for computing inharmonicity
    Output:
        meanInharm (float or np.float): mean inharmonicity over all the frames between the time interval 
                                        t1 and t2. 
    """
    # 0. Read the audio file and obtain an analysis window
    fs, x = UF.wavread(inputFile)
    w = get_window(window, M)
    # 1. Use harmonic model to compute the harmonic frequencies and magnitudes
    # xhfreq, xhmag, xhphase = harmonicModelAnal(x, fs, w, N, H, t, nH, minf0, maxf0, f0et, harmDevSlope=0.01, minSineDur=.02)
    xhfreq, xhmag, xhphase = HM.harmonicModelAnal(x, fs, w, N, H, t, nH, minf0, maxf0, f0et)

    # 2. Extract the time segment in which you need to compute the inharmonicity.
    l1 = int(round((t1 * fs + 1 / 2 * H) / H, 0))  # estemated frame ate t1
    l2 = int(round((t2 * fs + 1 / 2 * H) / H, 0))  # estemated frame ate t2
    xhfreq = xhfreq[l1:l2]
    xhmag = xhmag[l1:l2]
    xhphase = xhphase[l1:l2]

    # 3. Compute the mean inharmonicity of the segment
    r = np.arange(1, nH + 1)
    R = nH
    #r = np.tile(r,(xhfreq.size,1))
    # fr = r*f0 sqrt(1 + Br2)
    I = []
    for Ival in range(xhfreq.shape[0]):
        temp = (np.abs(xhfreq[Ival] - r * xhfreq[Ival, 0]))
        I = np.append(I, 1 / R * np.sum(temp / r))
    meanInharm = 1 / (l2 - l1 + 1) * np.sum(I)
    return(meanInharm)


if __name__ is "__main__":
    #I, meanInharm, xhfreq = estimateInharmonicity()
    #print(meanInharm)
    # Test case 1: If you run your code with inputFile = '../../sounds/piano.wav', t1=0.2, t2=0.4, window='hamming', M=2047, N=2048, H=128, f0et=5.0, t=-90, minf0=130, maxf0=180, nH = 25,
    # the returned output should be 1.4543.
    # meanInharm = estimateInharmonicity(inputFile = '../../sounds/piano.wav', t1=0.2, t2=0.4, window='hamming', M=2047, N=2048, H=128, f0et=5.0, t=-90, minf0=130, maxf0=180, nH=25)
    meanInharm = estimateInharmonicity(inputFile='../sms-tools/sounds/piano.wav', t1=0.2, t2=0.4, window='hamming', M=2047, N=2048, H=128, f0et=5.0, t=-90, minf0=130, maxf0=180, nH=25)
    print(meanInharm, 1.4543, meanInharm / 1.4543)

    # Test case 2: If you run your code with inputFile = '../../sounds/piano.wav', t1=2.3, t2=2.55, window='hamming', M=2047, N=2048, H=128, f0et=5.0, t=-90, minf0=230, maxf0=290, nH = 15,
    # the returned output should be 1.4874.
    # meanInharm = estimateInharmonicity(inputFile = '../../sounds/piano.wav', t1=2.3, t2=2.55, window='hamming', M=2047, N=2048, H=128, f0et=5.0, t=-90, minf0=230, maxf0=290, nH=15)
    meanInharm = estimateInharmonicity(inputFile='../sms-tools/sounds/piano.wav', t1=2.3, t2=2.55, window='hamming', M=2047, N=2048, H=128, f0et=5.0, t=-90, minf0=230, maxf0=290, nH=15)
    print(meanInharm, 1.4874, meanInharm / 1.4874)

    # Test case 3: If you run your code with inputFile = '../../sounds/piano.wav', t1=2.55, t2=2.8, window='hamming', M=2047, N=2048, H=128, f0et=5.0, t=-90, minf0=230, maxf0=290, nH = 5
    # the returned output should be 0.1748.
    # meanInharm = estimateInharmonicity(inputFile = '../../sounds/piano.wav', t1=2.55, t2=2.8, window='hamming', M=2047, N=2048, H=128, f0et=5.0, t=-90, minf0=230, maxf0=290, nH = 5)
    meanInharm = estimateInharmonicity(inputFile='../sms-tools/sounds/piano.wav', t1=2.55, t2=2.8, window='hamming', M=2047, N=2048, H=128, f0et=5.0, t=-90, minf0=230, maxf0=290, nH = 5)
    print(meanInharm, 0.1748, meanInharm / 0.1748)