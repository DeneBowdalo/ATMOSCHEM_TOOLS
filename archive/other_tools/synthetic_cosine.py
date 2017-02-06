import numpy as np
import csv
import glob
import datetime as datetime
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import logging as log
import nappy
import scipy.stats as stats
import lomb_phase
import modules
import pylab
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq, rfft

#set plotting area & background to white
fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(111)	
time = np.arange(0,40,1./24)

amplitude = 40
time_offset = 0
ave_species = 50
period = 365.25

print 'Amplitude = ', amplitude
print 'Time offset = ', time_offset
print 'Vertical_Average_Shift = ', ave_species
print 'Period = ', period
	
cosine_wave = ave_species + (amplitude*np.cos(2*np.pi*(time-time_offset)/period))

amplitude = 10
time_offset = 0
ave_species = 50
period = 182.625

print 'Amplitude = ', amplitude
print 'Time offset = ', time_offset
print 'Vertical_Average_Shift = ', ave_species
print 'Period = ', period

cosine_wave_2 = ave_species + (amplitude*np.sin(2*np.pi*(time-time_offset)/period))

amplitude = 20
time_offset = 0
ave_species = 50
period = 1

print 'Amplitude = ', amplitude
print 'Time offset = ', time_offset
print 'Vertical_Average_Shift = ', ave_species
print 'Period = ', period

cosine_wave_3 = ave_species + (amplitude*np.cos(2*np.pi*(time-time_offset)/period))

cosine_wave = cosine_wave_3#+cosine_wave_2+cosine_wave_3

# Remove annual period
	

#Remove half- annual signal


#Remove diurnal signal

#ofac = 2.00001
samp_freq = 24

#Normalise data
standard_deviation_model_p = np.std(cosine_wave)
mean_model_p = np.mean(cosine_wave)
normal_model = cosine_wave-mean_model_p
normal_model = normal_model/standard_deviation_model_p

#FFT 
samp_spacing=float(1./24.)
fft_model=np.abs(fft(cosine_wave))#**2 #calculate the magnitude of signal components
frequencies = fftfreq(cosine_wave.size, d = samp_spacing)
keep_freq = frequencies>0
frequencies, fft_model = frequencies[keep_freq], fft_model[keep_freq]
model_periods = 1 /frequencies
fft_model = fft_model/len(cosine_wave)
#fft_model = fft_model/samp_freq


#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
#nyquist_freq_lomb_model = frequencies[-1]
#Si_lomb_model = np.mean(fft_model)*nyquist_freq_lomb_model
#print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2

#print 'Time series amplitude = ', np.var(cosine_wave)

#Model lomb
fx, fy, vari, F = lomb_phase.lomb(time,cosine_wave)

#print F

magnitude = np.abs(F[1:])

#Divide output by sampling frequency
#fy = fy/samp_freq

#get amplitude
#amplitude = [np.sqrt(a**2 + b**2) for a,b in zip(real,imag)]
#amplitude = np.array(amplitude)
#amplitude = np.sqrt(fy*2*vari*2/float(len(cosine_wave)))
magnitude = magnitude/len(cosine_wave)


#print amplitude
#convert PSD to ESD
#fy = fy*len(cosine_wave)
#fy = np.sqrt(fy)


#lomb_amplitude = np.sum(amplitude)
#print lomb_amplitude
#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
#nyquist_freq_lomb_model = fx[-1]
#Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
#print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 

#Plot up
plt.plot(1./fx,magnitude,color='blue', label='Idealised Annual Cosine wave')	
plt.grid(True)
leg=plt.legend(loc=4, prop={'size':21})
#leg.get_frame().set_alpha(0.4)
plt.xlabel('Days', fontsize=21)
plt.ylabel('O3 Conc. (ppbV)', fontsize=21)
plt.title('Annual O3 Cosine wave ',fontsize=21)
#ax.text(0, 0.99,'Ratio of Actual Wave to Idealised Cosine Wave = %10.7f' %(np.average(difference)), ha='left', va='center', transform=ax.transAxes, fontweight ='bold')

#plt.ylim(-4,4)
#plt.xlim(0,10)
plt.show()

#set plotting area & background to white
#fig=plt.figure(figsize=(20,12))
#fig.patch.set_facecolor('white')
#ax = fig.add_subplot(111)

#annual_amplitude_test = fft_model
#fft_annual_amplitude = fft_model[annual_amplitude_test]

#plt.xlim(-1,1)

#frequencies = frequencies[1:]
#fft_model = fft_model[1:]
#fft_amplitude = np.sum(np.abs(fft_model))
#print 'FFT Amplitude = ', fft_amplitude
#print 'FFT Amplitude corrected = ', fft_amplitude*2

#plt.plot(frequencies, fft_model , color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='FFT')
#plt.xlabel('Frequency (1/Days)', fontsize=21)
#plt.ylabel('Magnitude (ppbV)', fontsize=21)
#plt.title('FFT Transformed Daily Signal ',fontsize=21)
#plt.plot(fx, fy , color = 'red', marker = 'x', alpha = 0.75,markersize=2, label='Lomb')

#plt.ylim(1e-4,800)

#plt.show()









