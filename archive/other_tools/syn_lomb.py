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

b = np.arange(0,550.5,hourly_step)

#b = b[:-1]

print len(b)

pi2 = np.pi*2


p1 = 1.
p2 = 90.
amp1 = 3.
amp2 = 10.
phase1 = 0.
phase2 = 0.

cos_waveform1 = amp1*(np.cos((pi2*b/p1)-(phase1)))
cos_waveform2 = amp2*(np.cos((pi2*b/p2)-(phase2)))

vals = 50 + (cos_waveform1+cos_waveform2)

print '\nActual Amplitude 1 = %s'%(amp1)
print 'Actual Phase 1 = %s'%(phase1)
print 'Actual Amplitude 2 = %s'%(amp2)
print 'Actual Phase 2 = %s'%(phase2)
print '------------------------------------------\n'

noise = np.random.normal(0,5,len(vals))
vals = vals+noise

b = np.array(b)
vals = np.array(vals)
orig_b = np.copy(b)
orig_vals = np.copy(vals)

plt.plot(b,vals)
plt.show()

#-------------------------------
#Put in gaps
all_frac_points = []
all_lomb_diff = []
all_fft_diff = []

for n in np.arange(100,99,-1.):
    print n
    
    n_points = (len(orig_vals)/100.)*n
    n_points = int(n_points)
    
    if n_points == 100:
        b = orig_b[:]
        vals = orig_vals[:]
    else:
        frac_points = random.sample(xrange(0,len(orig_vals)),n_points)
        frac_points = np.sort(frac_points)
        all_frac_points.append(n)
        b = orig_b[frac_points]
        vals = orig_vals[frac_points]
        
    #Lomb OUT OF BOX
    #-----------------------------------------------
    ofac = 1
    periods,mag,ph,fr,fi,amp_corr = modules.take_lomb(b,vals,ofac,hourly_step,w=False)
    
    closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
    closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
    daily_mag1 = mag[closest_period_1]
    daily_phase1 = ph[closest_period_1]
    daily_mag2 = mag[closest_period_2]
    daily_phase2 = ph[closest_period_2]

    print '\nLOMB OUT OF BOX'
    print 'Est. Amplitude 1 = ',daily_mag1
    print 'Est. Phase 1 = ',daily_phase1
    print 'Est. Amplitude 2 = ',daily_mag2
    print 'Est. Phase 2 = ',daily_phase2
    print '------------------------------------------\n'
       
    #Lomb WINDOWED
    #-----------------------------------------------
    ofac = 1
    periods,mag,ph,fr,fi,amp_corr = modules.take_lomb(b,vals,ofac,hourly_step,w=True)
    
    closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
    closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
    daily_mag1 = mag[closest_period_1]
    daily_phase1 = ph[closest_period_1]
    daily_mag2 = mag[closest_period_2]
    daily_phase2 = ph[closest_period_2]

    print '\nLOMB WINDOWED'
    print 'Est. Amplitude 1 = ',daily_mag1
    print 'Est. Phase 1 = ',daily_phase1
    print 'Est. Amplitude 2 = ',daily_mag2
    print 'Est. Phase 2 = ',daily_phase2
    print '------------------------------------------\n'
    
    #Lomb WINDOWED, OVERSAMPLED
    #-----------------------------------------------
    ofac = 4
    periods,mag,ph,fr,fi,amp_corr = modules.take_lomb(b,vals,ofac,hourly_step,w=True)
    
    closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
    closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
    daily_mag1 = mag[closest_period_1]
    daily_phase1 = ph[closest_period_1]
    daily_mag2 = mag[closest_period_2]
    daily_phase2 = ph[closest_period_2]

    print '\nLOMB WINDOWED, OVERSAMPLED'
    print 'Est. Amplitude 1 = ',daily_mag1
    print 'Est. Phase 1 = ',daily_phase1
    print 'Est. Amplitude 2 = ',daily_mag2
    print 'Est. Phase 2 = ',daily_phase2
    print '------------------------------------------\n'
    
    #Lomb WINDOWED, OVERSAMPLED, INTERPOLATED
    #-----------------------------------------------
    ofac = 4
    periods,mag,ph,fr,fi,amp_corr = modules.take_lomb(b,vals,ofac,hourly_step,w=True)   
    closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
    closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2)) 
    daily_mag1,daily_phase1 = modules.periodic_interp(fr,fi,1000,periods,p1,len(vals),amp_corr,window='hanning')
    daily_mag2,daily_phase2 = modules.periodic_interp(fr,fi,1000,periods,p2,len(vals),amp_corr,window='hanning')                                                                                                    
    print '\nLOMB WINDOWED, OVERSAMPLED, INTERPOLATED'
    print 'Est. Amplitude 1 = ',daily_mag1
    print 'Est. Phase 1 = ',daily_phase1
    print 'Est. Amplitude 2 = ',daily_mag2
    print 'Est. Phase 2 = ',daily_phase2
    print '------------------------------------------\n'
    
    #Lomb WINDOWED, SPECIFIC FREQUENCIES
    #-----------------------------------------------
    
    periods,mag,ph,fr,fi = modules.take_lomb_spec(b,vals,w=True,key_periods=[p1,p2])
    
    closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
    closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))

    daily_mag1 = mag[closest_period_1]
    daily_phase1 = ph[closest_period_1]
    daily_mag2 = mag[closest_period_2]
    daily_phase2 = ph[closest_period_2]

    print '\nLOMB WINDOWED, SPECIFIC FREQUENCIES'

    print 'Est. Amplitude 1 = ',daily_mag1
    print 'Est. Phase 1 = ',daily_phase1
    print 'Est. Amplitude 2 = ',daily_mag2
    print 'Est. Phase 2 = ',daily_phase2
    print '------------------------------------------\n'
                                                                                            

    #FFT - OUT OF BOX
    #----------------------------------------------------------------
    ofac=1
    periods,mag,ph,fr,fi,fft_array,amp_corr = modules.take_fft(b,vals,hourly_step,ofac,w=False)
    
    closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
    closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
    daily_mag1 = mag[closest_period_1]
    daily_phase1 = ph[closest_period_1]
    daily_mag2 = mag[closest_period_2]
    daily_phase2 = ph[closest_period_2] 
    
    print '\nFFT OUT OF BOX'
    print 'Est. Amplitude 1 = ',daily_mag1
    print 'Est. Phase 1 = ',daily_phase1
    print 'Est. Amplitude 2 = ',daily_mag2
    print 'Est. Phase 2 = ',daily_phase2
    print '------------------------------------------\n'
    
    #FFT - WINDOWED
    #----------------------------------------------------------------
    ofac=1
    periods,mag,ph,fr,fi,fft_array,amp_corr = modules.take_fft(b,vals,hourly_step,ofac,w=True)
    
    closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
    closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
    daily_mag1 = mag[closest_period_1]
    daily_phase1 = ph[closest_period_1]
    daily_mag2 = mag[closest_period_2]
    daily_phase2 = ph[closest_period_2] 
    
    print '\nFFT WINDOWED'
    print 'Est. Amplitude 1 = ',daily_mag1
    print 'Est. Phase 1 = ',daily_phase1
    print 'Est. Amplitude 2 = ',daily_mag2
    print 'Est. Phase 2 = ',daily_phase2
    print '------------------------------------------\n'
    
    #FFT - WINDOWED,OVERSAMPLED
    #----------------------------------------------------------------
    ofac=4
    periods,mag,ph,fr,fi,fft_array,amp_corr = modules.take_fft(b,vals,hourly_step,ofac,w=True)
    
    closest_period_1 = min(range(len(periods)), key=lambda i: abs(periods[i]-p1))
    closest_period_2 = min(range(len(periods)), key=lambda i: abs(periods[i]-p2))
    daily_mag1 = mag[closest_period_1]
    daily_phase1 = ph[closest_period_1]
    daily_mag2 = mag[closest_period_2]
    daily_phase2 = ph[closest_period_2] 
    
    print '\nFFT WINDOWED, OVERSAMPLED'
    print 'Est. Amplitude 1 = ',daily_mag1
    print 'Est. Phase 1 = ',daily_phase1
    print 'Est. Amplitude 2 = ',daily_mag2
    print 'Est. Phase 2 = ',daily_phase2
    print '------------------------------------------\n'

    #FFT - windowed,oversampled,interpolated
    #----------------------------------------------------------------
    #f = interpolate.interp1d(b, vals)
    #new_b = np.arange(np.min(b),np.max(b),hourly_step)
    #interp_vals =  f(new_b)
    
    ofac=4
    periods,mag,ph,fr,fi,fft_array,amp_corr = modules.take_fft(b,vals,hourly_step,ofac,w=True)
    
    daily_mag1,daily_phase1= modules.periodic_interp(fr,fi,1000,periods,p1,len(vals),amp_corr,window='hanning')
    daily_mag2,daily_phase2 = modules.periodic_interp(fr,fi,1000,periods,p2,len(vals),amp_corr,window='hanning')
    
    print '\nFFT WINDOWED, OVERSAMPLED, INTERPOLATED'
    print 'Est. Amplitude 1 = ',daily_mag1
    print 'Est. Phase 1 = ',daily_phase1
    print 'Est. Amplitude 2 = ',daily_mag2
    print 'Est. Phase 2 = ',daily_phase2
    print '------------------------------------------\n'
