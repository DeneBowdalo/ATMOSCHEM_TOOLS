import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
#from mpl_toolkits.basemap import Basemap
#from matplotlib.patches import Polygon
from cmath import *
from math import *

daily_period = 1
half_annual_period = (365.25/2)
annual_period = 365.25

f = 'MH_O3_SFC_03041987_29102013.npy'

#iterate through sites and take FFT

daily_mag_array = []
half_annual_mag_array = []
annual_mag_array = []

daily_phase_array = []
half_annual_phase_array = [] 
annual_phase_array = []

full_time = np.arange(0,9445,1./24)

read = np.load(f)
times = read[0,:]
read = read[1,:]
valid = read >= 0
read = read[valid]
times = times[valid]
	
#linearly interpolate values
read = np.interp(full_time,times,read)

#FFT 
samp_spacing=float(1./24.)
fft_calc=fft(read)#**2 #calculate the magnitude of signal components
fft_model = np.abs(fft_calc)
fft_phase = [phase(i) for i in fft_calc]
fft_phase = np.array(fft_phase)
model_frequencies = fftfreq(read.size, d = samp_spacing)
keep_freq = model_frequencies>0
model_frequencies, fft_model, fft_phase = model_frequencies[keep_freq], fft_model[keep_freq], fft_phase[keep_freq]
model_periods = 1 /model_frequencies

fft_model = fft_model/(len(read)/2)
fft_model = np.array(fft_model)
	
#calculations for mags and phases of key periods

closest_daily_period_index = min(range(len(model_periods)), key=lambda i: abs(model_periods[i]-daily_period))
closest_half_annual_period_index = min(range(len(model_periods)), key=lambda i: abs(model_periods[i]-half_annual_period))
closest_annual_period_index = min(range(len(model_periods)), key=lambda i: abs(model_periods[i]-annual_period))

daily_mag = fft_model[closest_daily_period_index]
half_annual_mag = fft_model[closest_half_annual_period_index]
annual_mag = fft_model[closest_annual_period_index]

daily_phase = fft_phase[closest_daily_period_index]
half_annual_phase = fft_phase[closest_half_annual_period_index]
annual_phase = fft_phase[closest_annual_period_index]

daily_mag_array = np.append(daily_mag_array,daily_mag)
half_annual_mag_array = np.append(half_annual_mag_array,half_annual_mag)
annual_mag_array = np.append(annual_mag_array,annual_mag)
	
daily_phase_array = np.append(daily_phase_array,daily_phase)
half_annual_phase_array = np.append(half_annual_phase_array,half_annual_phase)	
annual_phase_array = np.append(annual_phase_array,annual_phase)


fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')
	
plt.loglog(model_periods, fft_model, color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='Data')	
plt.show()

#save fft magitude
np.save('MH_obs_daily_magnitudes',daily_mag_array)
np.save('MH_obs_half_annual_magnitudes',half_annual_mag_array)
np.save('MH_obs_annual_magnitudes',annual_mag_array)

#save fft phase
np.save('MH_obs_daily_phases',daily_phase_array)	
np.save('MH_obs_half_annual_phases',half_annual_phase_array)
np.save('MH_obs_annual_phases',annual_phase_array)

