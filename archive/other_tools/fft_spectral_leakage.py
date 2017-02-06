import matplotlib.pyplot as plt 
from scipy.fftpack import fftfreq,fft
from scipy import signal
import numpy as np
from cmath import *
import operator
import modules
from multiprocessing import Pool
import time as now
from operator import itemgetter

#Time array in hourly intervals
hourly_step = 1./24

a = np.arange(0,365.25*6,hourly_step)

pi2 = np.pi*2

daily_amplitude = 2
half_annual_amplitude = 5
annual_amplitude = 10

cos_waveform1 = daily_amplitude*(np.cos((pi2*a/1)-(0)))
cos_waveform2 = half_annual_amplitude*(np.cos((pi2*a/(365.25/2))-0))
cos_waveform3 = annual_amplitude*(np.cos((pi2*a/365.25)+(np.pi/2)))

#waveform = (cos_waveform1)
waveform = 50 + (cos_waveform1 + cos_waveform2 + cos_waveform3)
#waveform = 50 + (cos_waveform1)

remove = np.arange(100,2100,1)
remove2 = np.arange(2500,4500,1)
remove3 = np.arange(6200,8200,1)
remove4 = np.arange(12000,14000,1)
remove5 = np.arange(19000,21000,1)
remove6 = np.arange(23000,25000,1)
remove7 = np.arange(30000,32000,1)

#for i in remove:
#	waveform[i] = -99

#for i in remove2:
#    waveform[i] = -99

#for i in remove3:
#    waveform[i] = -99

#for i in remove4:
#    waveform[i] = -99

#for i in remove5:
#    waveform[i] = -99

#for i in remove6:
#    waveform[i] = -99

#for i in remove7:
#    waveform[i] = -99
#first year/last year waveform

#do FFT
#test = waveform > 0
#waveform  = waveform[test]
#altered_a = a[test]

#waveform = np.interp(a,altered_a,waveform)
waveform_mean = np.mean(waveform)
waveform = waveform - waveform_mean

#window = signal.hamming(len(waveform))
window = np.kaiser(len(waveform),4)
waveform = waveform*window

#print len(waveform)

fft_array = [2,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288]

upper_i = np.searchsorted(fft_array,len(waveform))
gap =  fft_array[upper_i+1]-len(waveform)
append_array = [0]*gap

waveform=np.append(waveform,append_array)

fft_array = fft(waveform)
fft_cos_amp = fft_array.real
fft_sin_amp = fft_array.imag
fft_mag = np.abs(fft_array)
fft_phase = [phase(i) for i in fft_array]
fft_phase = np.array(fft_phase)
fft_freqs = fftfreq(len(waveform),hourly_step) 

valid = fft_freqs > 0
fft_mag, fft_phase, fft_freqs, fft_cos_amp, fft_sin_amp = fft_mag[valid], fft_phase[valid], fft_freqs[valid], fft_cos_amp[valid], fft_sin_amp[valid]

fft_periods = 1./fft_freqs
fft_cos_amp = fft_cos_amp/(len(a)/2)
fft_sin_amp = fft_sin_amp/(len(a)/2)
fft_mag = fft_mag/len(a)
fft_mag = fft_mag*2
#fft_mag = fft_mag* np.sqrt(8./3)
amp_corr = 1./(sum(window)/len(window))
#print 'amp corr', amp_corr
fft_mag = fft_mag * amp_corr

daily_index = np.argmin(np.abs(fft_periods - 1))
recon = fft_mag[daily_index]
#recon =  np.sqrt(fft_mag[daily_index]**2+fft_mag[daily_index-1]**2+fft_mag[daily_index+1]**2)
#recon =  np.sqrt(fft_mag[daily_index]**2+fft_mag[daily_index-1]**2+fft_mag[daily_index+1]**2+fft_mag[daily_index-2]**2+fft_mag[daily_index+2]**2)
#recon =  np.sqrt(fft_mag[daily_index]**2+fft_mag[daily_index-1]**2+fft_mag[daily_index+1]**2+fft_mag[daily_index-2]**2+fft_mag[daily_index+2]**2+fft_mag[daily_index-3]**2+fft_mag[daily_index+3]**2)
#recon =  np.sqrt(fft_mag[daily_index]**2+fft_mag[daily_index-1]**2+fft_mag[daily_index+1]**2+fft_mag[daily_index-2]**2+fft_mag[daily_index+2]**2+fft_mag[daily_index-3]**2+fft_mag[daily_index+3]**2+fft_mag[daily_index-4]**2+fft_mag[daily_index+4]**2)
#recon =  np.sqrt(fft_mag[daily_index]**2+fft_mag[daily_index-1]**2+fft_mag[daily_index+1]**2+fft_mag[daily_index-2]**2+fft_mag[daily_index+2]**2+fft_mag[daily_index-3]**2+fft_mag[daily_index+3]**2+fft_mag[daily_index-4]**2+fft_mag[daily_index+4]**2+fft_mag[daily_index-5]**2+fft_mag[daily_index+5]**2)
percent = (100/2) * recon 
print percent

#plot up
fig = plt.figure(figsize=(14, 10 ))
fig.patch.set_facecolor('white')
#ax1 = plt.subplot(211)
#plt.plot(a,waveform)
plt.semilogx(fft_periods, fft_mag)

plt.ylim(0,15)

#plt.xlim(0.99,1.01)

#plt.xlabel('Period (Days)', fontsize = 16)
#plt.ylabel('Magnitude (ppbV)', fontsize = 16)
#plt.title('FFT', fontsize = 16)

plt.grid()
#plt.yticks([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
plt.axvline(1, color = 'red', linestyle = '--')
plt.axvline(365.25/2, color = 'red', linestyle = '--')
plt.axvline(365.25, color = 'red', linestyle = '--')

plt.text(0.6,14,"Daily Period")
plt.text(50,13,"Half Annual Period")
plt.text(230,14,"Annual Period")

plt.show()
