import numpy as np
import csv
import glob
import datetime as datetime
import matplotlib.pyplot as plt
import logging as log
#import nappy
from cmath import *
import scipy.stats as stats
import lomb
import modules
from scipy.io import netcdf
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq, rfft, ifft
from pandas import *
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker



#select speci (or species) & location. Only Cape Verde has multiple species


location  = 'Cape_Verde'

if location == 'Cape_Verde':
    real_lat =-90 #16.848   #For CVO
    real_lon =-180 #-24.871

if location == 'Mauna_Loa':
	real_lat = 19.54
	real_lon = -155.58

species_list = ['O3']#,'CO','NO','NO2','C2H6','C3H8','DMS','TRA_6','ACET','GMAO_TEMP','GMAO_PSFC','GMAO_WIND','GMAO_RADSW','GMAO_RHUM']

	
def readfile(filename):
	read = np.load(filename)
	model_time = read[:,0]
	model = read[:]
	model_time = np.float64(model_time)
	model = read
	model_time = np.arange(len(model))
	model = np.float64(model)

	return model_time, model

def readfile2(filename):
    read = np.load(filename)
    model = read
    model_time = np.arange(len(model))
    model = np.float64(model)
    model = model *1e9
    model = model.flatten()	

    return model_time, model

def obs_O3_reader(filename):
    for files in filename:
        print files
        reader=csv.reader(open(files,'rb'), delimiter=' ', skipinitialspace=True)
        for row in reader:
            new = row[:]
            try:
                data.append(new)
            except:
                data=[new]

    data=np.array(data)
    year = data[:,1]
    month= data[:,2]
    day  = data[:,3]
    almost_date = [a+b for a,b in zip(year, month)]
    date = [a+b for a,b in zip(almost_date, day)]
    date = np.array(date)
    date = date.astype(int)
    hour = data[:,4]
    time = [i+'00' for i in hour]
    time = np.array(time)
    time = time.astype(int)
    vals = data[:,5]
    vals = np.float64(vals)
    return date, time, vals

#Get variables for location chosen
obsfile,loc_label,model_index = modules.location_check(location)

#Define sampling intervals
samp_spacing = 1./24.

counter = 0

for species in species_list:
	units, obs_data_name, unit_cut, species_type, actual_species_name, obs_switch, model_cut_switch, ofac = modules.obs_variable_finder(species)

#set plotting area & background to white
	fig=plt.figure(figsize=(20,12))
	fig.patch.set_facecolor('white')
	ax = plt.subplot(111)

	model , names = modules.readfile("binary_logs/GEOS_v90103_4x5_CV_logs.npy","001")

    # Processes the model date 
	date = model[:,0]
	time = model[:,1]
	model_time = modules.date_process(date,time)

	k=names.index('O3')
	model = model[:,k]*1e9



#Define sampling frequency
	samp_freq = 24

#FFT 
	samp_spacing=float(1./24.)

	#window = np.hamming(len(model))
	#wave_mean = np.mean(model)
	#model = model - wave_mean

	#model = model*window

	#fft_array = [2,8,16,32,64,128,256,512,1024,2048,4096,8192,16384,32768,65536,131072,262144,524288]

	#upper_i = np.searchsorted(fft_array,len(model))
	#gap =  fft_array[upper_i]-len(model)
	#append_array = [0]*gap

	#model = np.append(model,append_array)



	fft_array = fft(model)
	fft_model=np.abs(fft(model))#**2 #calculate the magnitude of signal components
	fft_real = fft_array.real
	fft_imag = fft_array.imag
	fft_phase = [phase(i) for i in fft_array]
	fft_phase = np.array(fft_phase)
	fft_freqs = fftfreq(model.size, d = samp_spacing)
	keep_freq = fft_freqs > 0
	cut_freqs,cut_FFT = fft_freqs[keep_freq], fft_model[keep_freq]
	model_periods = 1 /cut_freqs
	cut_FFT = cut_FFT/len(model)
	cut_FFT = cut_FFT*2
	#amp_corr = 1./(sum(model)/len(model))
	#fft_model = fft_model * amp_corr	
	#cut_FFT = cut_FFT * amp_corr
	#fft_model = fft_model/len(model)

	#output = ifft( fft( waveform ) )

	print 1./fft_freqs[1]

#Process fft array
	daily_real = np.copy(fft_real)
	daily_imag = np.copy(fft_imag)
	#half_annual_real = np.copy(fft_real)
	#half_annual_imag = np.copy(fft_imag)
	annual_real = np.copy(fft_real)
	annual_imag = np.copy(fft_imag)
	synoptic_noise_real = np.copy(fft_real)
	synoptic_noise_imag = np.copy(fft_imag)	
	weather_noise_real = np.copy(fft_real)
	weather_noise_imag = np.copy(fft_imag)
	medium_noise_real = np.copy(fft_real)
	medium_noise_imag = np.copy(fft_imag)
	long_noise_real = np.copy(fft_real)
	long_noise_imag = np.copy(fft_imag)

	pos_daily_index = np.argmin(np.abs(fft_freqs - (1.)))
	neg_daily_index = np.argmin(np.abs(fft_freqs - (-1.)))
	#pos_half_annual_index = np.argmin(np.abs(fft_freqs - (1./182.625)))
	#neg_half_annual_index = np.argmin(np.abs(fft_freqs - (1./-182.625)))
	pos_annual_index = np.argmin(np.abs(fft_freqs - (1./365)))
	neg_annual_index = np.argmin(np.abs(fft_freqs - (1./-365)))
	start_synoptic_noise_pos_index = np.argmin(np.abs(fft_freqs - (1./10))) 	
	end_synoptic_noise_pos_index = np.argmin(np.abs(fft_freqs - (1.)))	
	start_synoptic_noise_neg_index = np.argmin(np.abs(fft_freqs - (1./-10)))
	end_synoptic_noise_neg_index = np.argmin(np.abs(fft_freqs - (-1.)))	
	start_weather_noise_pos_index = pos_daily_index+1
	end_weather_noise_pos_index = np.argmin(np.abs(fft_freqs - (1./0.08333)))
	start_weather_noise_neg_index = neg_daily_index+1
 	end_weather_noise_neg_index = np.argmin(np.abs(fft_freqs - (1./-0.08333)))
	start_medium_noise_pos_index = pos_annual_index+1
	end_medium_noise_pos_index = np.argmin(np.abs(fft_freqs - (1./10)))
	start_medium_noise_neg_index = neg_annual_index+1
	end_medium_noise_neg_index = np.argmin(np.abs(fft_freqs - (1./-10)))
	start_long_noise_pos_index = np.argmin(np.abs(fft_freqs - (1./2191)))
	end_long_noise_pos_index = pos_annual_index
	start_long_noise_neg_index = np.argmin(np.abs(fft_freqs - (1./-2191)))
	end_long_noise_neg_index = neg_annual_index
	
	

	real_daily_indices = [pos_daily_index,neg_daily_index]
	imag_daily_indices=[pos_daily_index,neg_daily_index]
	#real_half_annual_indices = [pos_half_annual_index,neg_half_annual_index]
	#imag_half_annual_indices=[pos_half_annual_index,neg_half_annual_index]
	real_annual_indices = [pos_annual_index,neg_annual_index]
	imag_annual_indices=[pos_annual_index,neg_annual_index]
	real_synoptic_noise_indices = np.arange(start_synoptic_noise_pos_index,end_synoptic_noise_pos_index)
	real_synoptic_noise_indices = np.append(real_synoptic_noise_indices, np.arange(start_synoptic_noise_neg_index,end_synoptic_noise_neg_index))
	imag_synoptic_noise_indices = np.arange(start_synoptic_noise_pos_index,end_synoptic_noise_pos_index)
	imag_synoptic_noise_indices = np.append(imag_synoptic_noise_indices, np.arange(start_synoptic_noise_neg_index,end_synoptic_noise_neg_index))
	real_weather_noise_indices = np.arange(start_weather_noise_pos_index,end_weather_noise_pos_index)
	real_weather_noise_indices = np.append(real_weather_noise_indices, np.arange(start_weather_noise_neg_index,end_weather_noise_neg_index))
	imag_weather_noise_indices = np.arange(start_weather_noise_pos_index,end_weather_noise_pos_index)
	imag_weather_noise_indices = np.append(imag_weather_noise_indices, np.arange(start_weather_noise_neg_index,end_weather_noise_neg_index))
	real_medium_noise_indices = np.arange(start_medium_noise_pos_index,end_medium_noise_pos_index)
	real_medium_noise_indices = np.append(real_medium_noise_indices, np.arange(start_medium_noise_neg_index,end_medium_noise_neg_index))
	imag_medium_noise_indices = np.arange(start_medium_noise_pos_index,end_medium_noise_pos_index)
	imag_medium_noise_indices = np.append(imag_medium_noise_indices, np.arange(start_medium_noise_neg_index,end_medium_noise_neg_index))
	real_long_noise_indices = np.arange(start_long_noise_pos_index,end_long_noise_pos_index)
	real_long_noise_indices = np.append(real_long_noise_indices, np.arange(start_long_noise_neg_index,end_long_noise_neg_index))
	imag_long_noise_indices = np.arange(start_long_noise_pos_index,end_long_noise_pos_index)
	imag_long_noise_indices = np.append(imag_long_noise_indices, np.arange(start_long_noise_neg_index,end_long_noise_neg_index))


	
	
	mask = np.ones(fft_real.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	mask[real_daily_indices] = False
	daily_real[mask] = 0
	mask = np.ones(fft_imag.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	mask[imag_daily_indices] = False
	daily_imag[mask] = 0

	#mask = np.ones(fft_real.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	#mask[real_half_annual_indices] = False
	#half_annual_real[mask] = 0
	#mask = np.ones(fft_imag.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	#mask[imag_half_annual_indices] = False
	#half_annual_imag[mask] = 0
	
	mask = np.ones(fft_real.shape,dtype=bool) #np.ones_like(a,dtype=bool)	
	mask[real_annual_indices] = False
	annual_real[mask] = 0
	mask = np.ones(fft_imag.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	mask[imag_annual_indices] = False
	annual_imag[mask] = 0

	mask = np.ones(fft_real.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	mask[real_synoptic_noise_indices] = False
	synoptic_noise_real[mask] = 0
	mask = np.ones(fft_imag.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	mask[imag_synoptic_noise_indices] = False
	synoptic_noise_imag[mask] = 0

	mask = np.ones(fft_real.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	mask[real_weather_noise_indices] = False
	weather_noise_real[mask] = 0
	mask = np.ones(fft_imag.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	mask[imag_weather_noise_indices] = False
	weather_noise_imag[mask] = 0

	mask = np.ones(fft_real.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	mask[real_medium_noise_indices] = False
	medium_noise_real[mask] = 0
	mask = np.ones(fft_imag.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	mask[imag_medium_noise_indices] = False
	medium_noise_imag[mask] = 0

	mask = np.ones(fft_real.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	mask[real_long_noise_indices] = False
	long_noise_real[mask] = 0
	mask = np.ones(fft_imag.shape,dtype=bool) #np.ones_like(a,dtype=bool)
	mask[imag_long_noise_indices] = False
	long_noise_imag[mask] = 0





 	daily_array = np.vectorize(complex)(daily_real,daily_imag)
	#half_annual_array = np.vectorize(complex)(half_annual_real,half_annual_imag)
	annual_array = np.vectorize(complex)(annual_real,annual_imag)
	synoptic_noise_array = np.vectorize(complex)(synoptic_noise_real,synoptic_noise_imag)
	weather_noise_array = np.vectorize(complex)(weather_noise_real,weather_noise_imag)	
	medium_noise_array = np.vectorize(complex)(medium_noise_real,medium_noise_imag)
	long_noise_array = np.vectorize(complex)(long_noise_real,long_noise_imag)	
    


	daily_inverse =  ifft(daily_array)
	#half_annual_inverse =  ifft(half_annual_array)
	annual_inverse =  ifft(annual_array)
	synoptic_noise_inverse = ifft(synoptic_noise_array)
	weather_noise_inverse = ifft(weather_noise_array)	
	medium_noise_inverse = ifft(medium_noise_array)	
	long_noise_inverse = ifft(long_noise_array)

	dc_component =  fft_model[0]/len(model)

	#daily_inverse = daily_inverse - ((fft_model[pos_daily_index]/len(model))*2)	
	#full_time_series = dc_component + (daily_inverse + annual_inverse + synoptic_noise_inverse + weather_noise_inverse + medium_noise_inverse + long_noise_inverse)
	full_time_series = dc_component + (long_noise_inverse)
	#full_time_series  = dc_component + (long_noise_inverse) 

	inverse_time = np.arange(len(daily_inverse))
	inverse_time = inverse_time/24.

#split fft into different periods
	weather_test = model_periods < 1
	synoptic_test = (model_periods > 1) & (model_periods < 10) 
	mid_freq_test = (model_periods > 10) & (model_periods < 360)
	low_freq_test = model_periods > 370

	annual_index = np.argmin(np.abs(model_periods - (365)))
	daily_index = np.argmin(np.abs(model_periods - (1)))

	annual_indices = [annual_index-1,annual_index,annual_index+1]
	daily_indices = [daily_index-1,daily_index,daily_index+1]


	weather_periods = model_periods[weather_test]
	weather_FFT = cut_FFT[weather_test]
	synoptic_periods = model_periods[synoptic_test]
	synoptic_FFT = cut_FFT[synoptic_test]
	mid_freq_periods = model_periods[mid_freq_test]
	mid_freq_FFT = cut_FFT[mid_freq_test]
	low_freq_periods = model_periods[low_freq_test]
	low_freq_FFT = cut_FFT[low_freq_test]
	annual_periods = model_periods[annual_indices] 
	annual_FFT = cut_FFT[annual_indices]
	daily_periods = model_periods[daily_indices]
	daily_FFT = cut_FFT[daily_indices]

	

	#plt.loglog(weather_periods,weather_FFT,color='blue')
	#plt.loglog(synoptic_periods,synoptic_FFT,color='red')	
	#plt.loglog(mid_freq_periods,mid_freq_FFT,color='purple')
	#plt.loglog(low_freq_periods,low_freq_FFT,color='green')
	#plt.loglog(annual_periods,annual_FFT, color = 'orange')
	#plt.loglog(daily_periods,daily_FFT,  color='black')

	#print annual_periods
	#print annual_FFT

#plot up
	#plt.loglog(1./fx, fy, color = 'red', marker = 'x', alpha = 0.75,markersize=2, label='GFDL Historical O3')
	#plt.loglog(model_periods, fft_model, color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='Mauna Loa GFDL Historical O3')
	#ax1 = plt.subplot(211)
	#ax1.loglog(model_periods, (cut_FFT/len(model))*2)
	#plt.plot(model_time,model, color='red')
	#ax2 = plt.subplot(212)
	plt.plot(inverse_time,full_time_series,color='red')
	#plt.plot(inverse_time,model,color='red')
	#plt.loglog(model_periods,cut_FFT,color='purple')
	#plt.plot(model_time,model)
	#ave_time = np.arange(24)
	#plt.plot(ave_time,ave_model)
#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
#nyquist_freq_lomb_model = frequencies[-1]
#Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
#print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 

	#Reverse array for smoothing
	#reversed_fft_obs = fft_obs[::-1]
	#reversed_obs_periods = obs_periods[::-1]

	#cut_obs_periods,smoothed_obs = modules.ema_keep_peaks(reversed_fft_obs, reversed_obs_periods, 20, 1000,0.01,3,6)
	#smoothed_obs = np.exp(smoothed_obs)
	
	#plt.loglog(cut_obs_periods, smoothed_obs, color = 'red', marker = 'x', alpha = 0.75,markersize=2, label='Mauna Loa Obs. O3')

	def form2(x, pos):
		""" This function returns a string with 3 decimal places, given the input x"""
		return '%.2f' % x

	def form6(x, pos):
		""" This function returns a string with 3 decimal places, given the input x"""
		return '%.6f' % x


	xformatter = FuncFormatter(form2)
	yformatter = FuncFormatter(form6)
 

	#ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
	#plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))

	#ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
	#plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))


	plt.tick_params(axis='x', which='major', labelsize=14)
	plt.tick_params(axis='y', which='major', labelsize=14)
	plt.tick_params(direction='out', pad=0.00005)

	#plt.xlim(0,5)
	#plt.ylim(5,50)
	plt.grid(True)
	leg=plt.legend(loc=4, prop={'size':21})
	#leg.get_frame().set_alpha(0.4)
	plt.xlabel('Period (Days)', fontsize = 16, labelpad = 8)
	plt.ylabel('Magnitude (ppbV)', fontsize = 16, labelpad = 8)
	#plt.ylabel('Conc. (ppbV)', fontsize = 16, labelpad = 8)
	#title = plt.title(r'Low-Frequency Noise Time Series for Surface $O_3$ from GEOS-Chem 4x5 Model (2006-2012), (constructed using IFFT)',fontsize=22)
	#title = plt.title(r'Reconstructed Time Series for Surface $O_3$ from GEOS-Chem 4x5 Model (2006-2012), (constructed using IFFT)',fontsize=22)	
	title = plt.title(r'FFT Spectrum for Surface $O_3$ from GEOS-Chem 4x5 Model (2006-2012)',fontsize=22)
	#title = plt.title(r'Time Series of Surface $O_3$ from GEOS-Chem 4x5 Model at Cape Verde',fontsize=22)

	title.set_y(1.01)
	counter+=1

	#plt.savefig('plots/v90103_nestedeurope_1year_GEOS5/normal/%s_%s_%s_%s.png'%(actual_species_name,location,mversion,res),dpi=fig.dpi)

plt.show()










