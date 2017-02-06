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
import lomb_phase

daily_period = 1
half_annual_period = (365.25/2)
annual_period = 365.25

files = glob.glob('chunked_binary_files/MH*')
files.sort()

print files

start_year =1988
end_year = 1993

counter = -1
#iterate through sites and take FFT
for f in files:
	counter+=1
	if counter==0:
		end_days=1827
	else:
		end_days=1826		
	if counter == 3:
		counter=-1
		
	daily_mag_array = []
	half_annual_mag_array = []
	annual_mag_array = []

	daily_phase_array = []
	half_annual_phase_array = [] 
	annual_phase_array = []

	full_time = np.arange(0,end_days,1./24)

	read = np.load(f)
	times = read[0,:]
	read = read[1,:]
	valid = read >= 0
	read = read[valid]
	times = times[valid]
	
	#linearly interpolate values
	#read = np.interp(full_time,times,read)

	#FFT 
	#samp_spacing=float(1./24.)
	#fft_calc=fft(read)#**2 #calculate the magnitude of signal components
	#fft_model = np.abs(fft_calc)
	#fft_phase = [phase(i) for i in fft_calc]
	#fft_phase = np.array(fft_phase)
	#model_frequencies = fftfreq(read.size, d = samp_spacing)
	#keep_freq = model_frequencies>0
	#model_frequencies, fft_model, fft_phase = model_frequencies[keep_freq], fft_model[keep_freq], fft_phase[keep_freq]
	#model_periods = 1 /model_frequencies

	#fft_model = fft_model/(len(read)/2)
	#fft_model = np.array(fft_model)

	#LSP
	window = np.kaiser(len(read),4)
	obs_mean = np.mean(read)	
	read = read - obs_mean
	read = read*window

	NOUT = 0.5*1*1*len(read)
	NOUT = int(NOUT)

	fa, fb, fft_model, fft_phase= lomb_phase.lomb(times,read,NOUT)
	model_periods = 1./fa
	amp_corr = 1./(sum(window)/len(window))
	fft_model = fft_model * amp_corr

	
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

	print daily_mag	
	plt.loglog(model_periods, fft_model, color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='Data')	
	plt.show()

	#save fft magitude
	np.save('chunked_mags_phases/LSP_WINDOW/%s_%s/MH_obs_daily_magnitudes'%(start_year,end_year),daily_mag_array)
	np.save('chunked_mags_phases/LSP_WINDOW/%s_%s/MH_obs_half_annual_magnitudes'%(start_year,end_year),half_annual_mag_array)
	np.save('chunked_mags_phases/LSP_WINDOW/%s_%s/MH_obs_annual_magnitudes'%(start_year,end_year),annual_mag_array)

	#save fft phase
	np.save('chunked_mags_phases/LSP_WINDOW/%s_%s/MH_obs_daily_phases'%(start_year,end_year),daily_phase_array)	
	np.save('chunked_mags_phases/LSP_WINDOW/%s_%s/MH_obs_half_annual_phases'%(start_year,end_year),half_annual_phase_array)
	np.save('chunked_mags_phases/LSP_WINDOW/%s_%s/MH_obs_annual_phases'%(start_year,end_year),annual_phase_array)

	start_year+=1
	end_year+=1

