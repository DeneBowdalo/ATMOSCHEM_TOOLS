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



daily_period = 1

hourly_step = 1./24

a = np.arange(0,25,hourly_step)
b = np.arange(0,200,hourly_step)


pi2 = np.pi*2


daily_amplitude = 2
half_annual_amplitude = 5
annual_amplitude = 10

cos_waveform1 = daily_amplitude*(np.cos((pi2*b/1)-(0)))
cos_waveform2 = half_annual_amplitude*(np.cos((pi2*b/(365.25/2))-0))
cos_waveform3 = annual_amplitude*(np.sin((pi2*b/365.25)+(0)))

#cos_waveform1 = daily_amplitude*(np.cos((pi2*b/10)+(0)))
#cos_waveform1 = daily_amplitude*(np.cos((pi2*b/1)-(0)))
vals = 50 + (cos_waveform1+cos_waveform2+cos_waveform3)

window = np.hanning(len(vals))
mean = np.mean(vals)
vals = vals - mean
vals = vals*window


NOUT = 0.5*4*1*len(vals)
NOUT = int(NOUT)

fa, fb, mag, ph= lomb_phase.lomb(b,vals,NOUT)
periods = 1./fa
amp_corr = 1./(sum(window)/len(window))
mag = mag * amp_corr
	

print np.min(ph)
print np.max(ph)

closest_daily_period_index = min(range(len(periods)), key=lambda i: abs(periods[i]-daily_period))
    
print periods[closest_daily_period_index+1]
      
daily_mag = mag[closest_daily_period_index]
daily_phase = ph[closest_daily_period_index]

plt.plot(b,vals)
#plt.show()

print 'daily magnitude', daily_mag
print 'daily phase = ', daily_phase