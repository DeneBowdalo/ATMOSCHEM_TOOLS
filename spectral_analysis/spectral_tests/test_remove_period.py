import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from cmath import *
from math import *
import csv
import datetime
import lomb_phase
from scipy import signal
import multiprocessing
import datetime
import time
import modules
from scipy import interpolate
import random

hourly_step = 1./24

b = np.arange(0,214.5,hourly_step)

#b = b[:-1]

print len(b)

pi2 = np.pi*2


p1 = 1.
p2 = 50.
amp1 = 3.
amp2 = 10.
phase_1 = 0.
phase_2 = 0.

a1_ratio = 100./amp1 
a2_ratio = 100./amp2

cos_waveform1 = amp1*(np.cos((pi2*b/p1)-(phase_1)))
cos_waveform2 = amp2*(np.cos((pi2*b/p2)-(phase_2)))

vals = 50 + (cos_waveform1+cos_waveform2)

print '\nActual Amplitude 1 = %s'%(amp1)
print 'Actual Phase 1 = %s'%(phase_1)
print 'Actual Amplitude 2 = %s'%(amp2)
print 'Actual Phase 2 = %s'%(phase_2)
print '------------------------------------------\n'

noise = np.random.normal(0,2,len(vals))
vals = vals+noise

b = np.array(b)
vals = np.array(vals)

new_b,new_vals = modules.lomb_remove_periodicity(b,vals,w = True,key_periods=[1.,50.])

plt.plot(b,vals)
plt.plot(new_b,new_vals)
plt.show()