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
from scipy import signal
import multiprocessing
import datetime
import time
import modules
import redfit
import random
import numpy.fft
from cmath import *
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
import scipy
import numpy.fft as FFT

daily_period = 1
annual_period = 365.25

hourly_step = 1./24

b = np.arange(200.75,500,hourly_step)

pi2 = np.pi*2


daily_amplitude = 10
half_annual_amplitude = 1
annual_amplitude = 1


#shift times to start from 0
b = modules.phase_start_correct(b)
print 'First time = ',b[0]
#b = b-b[0]
print b

cos_waveform1 = daily_amplitude*(np.cos((pi2*b/1.1)-(np.pi/2)))
#cos_waveform2 = half_annual_amplitude*(np.cos((pi2*b/(365.25/2))-0))
cos_waveform3 = annual_amplitude*(np.cos((pi2*b/365.25)-(0)))
#cos_waveform4 = daily_amplitude*(np.cos((pi2*b/1.0001)-(np.pi)))

#vals = 50 + (cos_waveform1+cos_waveform2+cos_waveform3)

vals = 50 + (cos_waveform1+cos_waveform3)

total_len = len(vals)

#gap data
remove_ints = []
for x in range(500):
    remove_ints.append(random.randint(0,len(vals)))
    remove_ints.append(x)

#print remove_ints
vals = np.delete(vals,remove_ints)
b = np.delete(b,remove_ints)

print '%s percent of dataset complete'%((100./total_len)*len(vals))

NOUT = 0.5*4*1*len(vals)
NOUT = int(NOUT)

periods, mag, ph,fr,fi,amp_corr= modules.take_lomb(b,vals,4,hourly_step)

closest_period_index = min(range(len(periods)), key=lambda i: abs(periods[i]-1.1))

print 'closest period = ', periods[closest_period_index]
period_closest = periods[closest_period_index]

#peak_raw = fft_raw[closest_period_index-10:closest_period_index+11]

peak_fr = fr[closest_period_index-5:closest_period_index+6]
peak_fi = fi[closest_period_index-5:closest_period_index+6]
peak_periods = periods[closest_period_index-5:closest_period_index+6]
peak_mag = mag[closest_period_index-5:closest_period_index+6]
peak_phase = ph[closest_period_index-5:closest_period_index+6]

closest_period_index_2 = min(range(len(peak_periods)), key=lambda i: abs(peak_periods[i]-1.1))


print 'peak periods =', peak_periods

peak_real = peak_fr/len(vals)
peak_real = peak_real*2
peak_real = peak_real*amp_corr

peak_imag = peak_fi/len(vals)
peak_imag = peak_imag*2
peak_imag = peak_imag*amp_corr

print 'Initial real = ',peak_real[closest_period_index_2]
print 'Initial imag = ',peak_imag[closest_period_index_2]
print 'Initial Amp = ', peak_mag[closest_period_index_2]
print 'Initial Phase = ', peak_phase[closest_period_index_2]

reverse_peak_periods = peak_periods[::-1]
reverse_peak_fr = peak_real[::-1]
reverse_peak_fi = peak_imag[::-1]

plt.plot(reverse_peak_periods,reverse_peak_fi,linestyle='None',marker='x',markersize = 14)

f = interp1d(reverse_peak_periods, reverse_peak_fr,kind='cubic')
f2 = interp1d(reverse_peak_periods, reverse_peak_fi,kind='cubic')
new_periods = np.linspace(np.min(reverse_peak_periods), np.max(reverse_peak_periods), 1000)
interp_fr = f(new_periods)
interp_fi = f2(new_periods)

closest_period_index_3 = min(range(len(new_periods)), key=lambda i: abs(new_periods[i]-1.1))

interp_fr_peak = interp_fr[closest_period_index_3] 
interp_fi_peak = interp_fi[closest_period_index_3]

plt.plot(new_periods,interp_fi,label='Cubic Spline Interpolation')

print 'Cubic Spline Interpolated fr = ',interp_fr_peak
print 'Cubic Spline Interpolated fi = ',interp_fi_peak

print 'Cubic Spline Interpolated Amp = ',np.abs(complex(interp_fr_peak,interp_fi_peak))
print 'Cubic Spline Interpolated Phase = ',np.angle(complex(interp_fr_peak,interp_fi_peak))


zoomfact = 998

interp_amp,interp_phase,newreal,newx = modules.periodic_interp(fr,fi,zoomfact,periods,1.1,len(vals),amp_corr)

plt.axvline(1.1)
plt.plot(newx,newreal,label='Hann Window Interpolation')
plt.legend()

plt.show()

print 'Hann Interpolated Amp = ',interp_amp
print 'Hann Interpolated Phase = ',interp_phase



#print 'yeah =', np.abs(complex(peak_real_max,peak_imag_max))

#peak_real = peak_fr
#peak_imag = peak_fi

#peak_real_correct = parabolic(peak_real,np.argmax(peak_real))
#peak_imag_correct = parabolic(peak_imag,10)

#peak_real_correct = peak_real_correct[1]
#peak_imag_correct = peak_imag_correct[1]

#print peak_real_correct
#print peak_imag_correct

#peak_complex = complex(peak_real_correct,peak_imag_correct)
#peak_correct_mag = np.abs(peak_complex)
#peak_correct_ph = np.angle(peak_complex)

#print 'corrected mag = ', peak_correct_mag
#print 'corrected ph = ', peak_correct_ph


#peak_mag = np.abs(peak_raw)
#peak_mag = (peak_mag/len(vals))*2.
#peak_phase = np.angle(peak_raw)

#var_mean = np.average(vals)



#Section uses ifft to generate time series of detrended, key period treated data
# make complex Fourier spectrum which corresponds to the Lomb-Scargle periodogram:
#F = [0]*((len(fr)*2)+1)
#set first real value to average 
#F[0] = complex(var_mean,0)

#Get reverse real and imaginary values
#rev_fr=fr[::-1] 
#rev_fi=fi[::-1]

#Fill Fourier Spectrum real and imaginary values
#for i in range(len(fr)):
#    F[i+1] = complex(fr[i],fi[i])
#for i in range(len(fr),len(fr)*2):
#    F[i+1] = complex(fr[i-len(fr)],-fi[i-len(fr)])
#F = np.array(F)

#Take ifft and just take real values
#output_ts = numpy.fft.ifft(F)
#output_ts = output_ts.astype('float64')
#time = np.arange(len(output_ts))
#output_ts = output_ts/np.hanning(len(output_ts))
#time = b[-len(vals):]
#output_ts = output_ts[-len(vals):]
#time = time[350:-350]
#output_ts = output_ts[350:-350]
#output_ts = output_ts+var_mean


#closest_period_index = min(range(len(periods)), key=lambda i: abs(periods[i]-1))

#closest_mag = mag[closest_period_index]
#up_mag = mag[closest_period_index-1]
#down_mag = mag[closest_period_index+1]

#closest_period = periods[closest_period_index]
#closest_ph = ph[closest_period_index]
#closest_fr = fr[closest_period_index]
#up_fr = fr[closest_period_index-1]
#down_fr = fr[closest_period_index+1]
#closest_fi = fi[closest_period_index]
#up_fi = fi[closest_period_index-1]
#down_fi = fi[closest_period_index+1]

#F = complex(closest_fr,closest_fi)

#window = signal.hanning(len(vals))
#amp_corr = 1./(sum(window)/len(window))


#recon_mag = np.abs(F)
#recon_mag = recon_mag/len(b)
#recon_mag = recon_mag*2
#recon_mag = recon_mag * amp_corr

#recon_ph = np.arctan2(closest_fi,closest_fr)

#print 'closest period = ',closest_period
#print 'all fr =', fr[closest_period_index-5:closest_period_index+5]
#print 'all fi =', fi[closest_period_index-5:closest_period_index+5]
#print 'Closest fr = ',closest_fr
#print 'Closest fi = ',closest_fi
#print 'Fr Up = ', up_fr
#print 'Fr Down = ', down_fr
#print 'Fi Up = ', up_fi
#print 'Fi Down = ', down_fi


#print 'Initial Mag = ', closest_mag 
#print 'Mag Up = ', up_mag
#print 'Mag Down = ', down_mag
#print 'Reconstructed Mag = ', recon_mag

#print 'Reconstructed Phase = ', recon_ph

#new_mag=[up_mag,recon_mag,down_mag]
#new_mag = modules.parabolic_amplitude_correct(mag,1., periods)


#correct_phase = modules.phase_offset_correct(1,periods,ph)
#correct_phase = modules.phase_start_point_correct(1.,correct_phase,b)
#print 'correct phase = ', correct_phase

#correct_phase = 3.404

#correct_phase = modules.convert_phase_units_actual_single(correct_phase,24)
#print 'correct phase = ', correct_phase
#annual_phase = modules.phase_offset_correct(1.,periods,ph)

#daily_phase = modules.phase_start_point_correct(1.,daily_phase,b)
#annual_phase = modules.phase_start_point_correct(365.25,annual_phase,b)

#print 'daily_phase = ', daily_phase	
#print 'annual_phase = ', annual_phase

#cos_waveform3 = 50+(daily_mag*(np.cos((pi2*b/365.25)-(annual_phase))))

#fig=plt.figure(figsize=(20,12))

#fig=plt.subplot(411)
#plt.plot(peak_periods,peak_real,linestyle='--',marker='x')

#fig=plt.subplot(412)
#plt.plot(peak_periods,peak_imag,linestyle='--',marker='x')

#fig=plt.subplot(413)
#plt.plot(peak_periods,peak_mag,linestyle='--',marker='x')

#plt.plot(period_closest,peak_mag,'x')

#fig=plt.subplot(414)
#plt.plot(peak_periods,peak_phase,linestyle='--',marker='x')


#plt.semilogx(periods,ph,color='black')
#plt.xlim(0,5)
#plt.ylim(48,52)
#plt.loglog(fft_periods,fft_mag,color='red')
#plt.plot(b,vals)
#plt.plot(b,cos_waveform3)
#plt.plot(time,output_ts)
#plt.show()
