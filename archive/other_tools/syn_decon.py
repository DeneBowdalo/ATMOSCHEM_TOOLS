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

b = np.arange(0,410,hourly_step)

#b = b[:-1]

pi2 = np.pi*2

p1 = 30.
p2 = 50.
amp1 = 21.
amp2 = 10.
phase1 = 0.
phase2 = 0.

cos_waveform1 = amp1*(np.cos((pi2*b/p1)-(phase1)))
cos_waveform2 = amp2*(np.cos((pi2*b/p2)-(phase2)))

vals = 50 + (cos_waveform1+cos_waveform2)

#plt.plot(b,vals,marker='x')
#noise = np.random.normal(0,5,len(vals))
#vals = vals+noise

orig_b = np.copy(b)

n_points = (len(vals)/100.)*95.
n_points = int(n_points)    
frac_points = random.sample(xrange(0,len(vals)),n_points)
frac_points = np.sort(frac_points)
non_frac_points = np.delete(np.arange(len(orig_b)),frac_points)
b = b[frac_points]
vals = vals[frac_points]

plt.plot(b,vals,marker='x')

b = np.array(b)
vals = np.array(vals)


#Lomb, out of box
#-----------------------------------------------
ofac = 4
#periods,mag,ph,fr,fi,amp_corr = modules.take_lomb(b,vals,ofac,hourly_step,w=False,kp=[p1,p2]) 
periods,mag,ph,fr,fi = modules.take_lomb_spec(b,vals,w=True,key_periods=[p1,p2])

closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
daily_mag1 = mag[closest_period_1]
daily_phase1 = ph[closest_period_1]
daily_mag2 = mag[closest_period_2]
daily_phase2 = ph[closest_period_2] 

print '\nLomb'
print 'Est. Amplitude 1 = ',daily_mag1
print 'Est. Phase 1 = ',daily_phase1
print 'Est. Amplitude 2 = ',daily_mag2
print 'Est. Phase 2 = ',daily_phase2
print '------------------------------------------\n'

ifft_time, ifft_ts = modules.lomb_ifft_spec(b,vals,periods,mag,ph,hourly_step)

#ifft_time,ifft_ts = modules.lomb_ifft([p1,p2],b,vals,periods,fr,fi,ofac,amp_corr,w=False)



plt.plot(ifft_time,ifft_ts)

#ftt
#------------------------------
ofac=4
periods,mag,ph,fr,fi,fft_array,amp_corr = modules.take_fft(b,vals,hourly_step,ofac,w=True)

closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
daily_mag1 = mag[closest_period_1]
daily_phase1 = ph[closest_period_1]
daily_mag2 = mag[closest_period_2]
daily_phase2 = ph[closest_period_2] 

print '\nFFT'
print 'Est. Amplitude 1 = ',daily_mag1
print 'Est. Phase 1 = ',daily_phase1
print 'Est. Amplitude 2 = ',daily_mag2
print 'Est. Phase 2 = ',daily_phase2
print '------------------------------------------\n'
ifft_time,ifft_ts = modules.fft_ifft(b,vals,fr,fi,fft_array,ofac,amp_corr,w=True)

#plt.plot(ifft_time,ifft_ts)
plt.show()















