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
import lomb_phase
import lomb_phase_spec

#TAKE SPECIFIC LOMB-SCARGLE - NEED FORTRAN FILE COMPILED WITH F2PY (lomb_phase_spec.so)
def take_lomb_spec(time,var,w = False,key_periods=[1./4.,1./3.,1./2.,1.,365.25/4.,365.25/3.,365.25/2.,365.25]):
    if w == True:
        window = signal.hanning(len(var))
        var_mean = np.mean(var)
        var = var - var_mean
        var = var*window
        amp_corr = 1./(sum(window)/len(window))
	
    freqs = 1./np.array(key_periods)
    
	#take lomb
    fb, mag, ph, fr, fi = lomb_phase_spec.lomb(time,var,freqs)
    
    
    if w == True:
        mag = mag * amp_corr
        
    periods = 1./freqs
    
    #CORRECT IMAGINARY COMPONENTS TO FFT EQUIVALENTS
    for i in range(len(fi)):
        if fi[i] < 0:
            fi[i] = fi[i]*-1
        elif fi[i] > 0:
            fi[i] = -fi[i]
    	
    return periods,mag,ph,fr,fi
    
#SUPERPOSE FUNDAMENTAL AND HARMONIC WAVEFORMS
def period_convolution(key_periods,time,mag,ph,mean):
    pi2 = np.pi*2
    
    if np.max(key_periods) == 1.:
        max_time = 1.
    elif np.max(key_periods) == 365.25:
        max_time = 365.25
    else:
        max_time = np.max(key_periods)
    
    for i in range(len(key_periods)):
        try:
            waveform = waveform + (mag[i]*(np.cos((pi2*time/key_periods[i])-(ph[i]))))
        except:
            waveform = mag[i]*(np.cos((pi2*time/key_periods[i])-(ph[i])))
    
    ff = np.sqrt(np.mean(np.square(waveform)))/np.mean(np.abs(waveform))
    waveform = mean+waveform
    waveform_min = np.min(waveform)
    waveform_max = np.max(waveform)
    waveform_min_ind = np.argmin(waveform)
    waveform_max_ind = np.argmax(waveform)
    mag = (waveform_max-waveform_min)/2.
    ph_min = time[waveform_min_ind]
    ph_max = time[waveform_max_ind]
    ratio = (pi2)/max_time
    ph_min = ratio*ph_min
    ph_max = ratio*ph_max
    
    return mag,ph_min,ph_max,waveform,ff

#TAKE LOMB AT SPECIFIC A PRIORI KNOWN FREQUENCIES AND SUBTRACT SUPERPOSED
#WAVEFORM FROM RAW TIME SERIES
def lomb_remove_periodicity(time,var,w = True,key_periods=[365.25/4.,365.25/3.,365.25/2.,365.25]):
    valid = var > 0
    cut_time = time[valid]
    cut_var = var[valid]
    mean = np.average(cut_var)
    
    periods,mag,ph,fr,fi = take_lomb_spec(cut_time,cut_var,w,key_periods)
    
    full_mag,full_min_ph,full_max_ph,full_waveform,full_ff = period_convolution(key_periods,time,mag,ph,mean)
    
    var[~valid] = np.NaN
    var = (var - full_waveform)+mean
    return time,var,full_waveform

#--------------------------------------------------------
#--------------------------------------------------------
#create synthetic data
hourly_step = 1./24
pi2 = np.pi*2
time = np.arange(0,214.5,hourly_step)

p1 = 1.
p2 = 50.
amp1 = 3.
amp2 = 10.
phase_1 = 0.
phase_2 = 0.

a1_ratio = 100./amp1 
a2_ratio = 100./amp2

cos_waveform1 = amp1*(np.cos((pi2*time/p1)-(phase_1)))
cos_waveform2 = amp2*(np.cos((pi2*time/p2)-(phase_2)))

#add average 
vals = 50 + (cos_waveform1+cos_waveform2)
#add noise
noise = np.random.normal(0,4,len(vals))
vals = vals+noise

time,new_vals,removed_waveform = lomb_remove_periodicity(time,vals,key_periods=[1.,50.])
plt.plot(time,vals,label='With Periodicity')
plt.plot(time,new_vals,label='Without Periodicity')
plt.plot(time,removed_waveform,label='Removed Periodicity')
plt.legend()
plt.show()