#fft_sin_amp = fft_sin_amp/(len(a)/2)
import numpy as np
import csv
import glob
import datetime as datetime
import matplotlib.pyplot as plt
import logging as log
from cmath import *
from scipy import signal
import scipy.stats as stats
import lomb
import lomb_phase
import modules
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq, rfft
import scipy.signal.spectral as spectral
from multiprocessing import Pool
import random
import time as now
import matplotlib.ticker
from matplotlib.ticker import FuncFormatter
from numba.decorators import jit, autojit

#select speci (or species) & location. Only Cape Verde has multiple species

species_list = ['O3']#,'CO','NO','NO2','C2H6','C3H8','DMS','TRA_6','ACET','GMAO_TEMP','GMAO_PSFC','GMAO_WIND','GMAO_RADSW','GMAO_RHUM']

#Location list for O3
#Arrival_Heights
#Barrow
#Lauder
#Mace_Head
#Mauna_Loa'
#Niwot_Ridge
#Ragged_Point
#South_Pole
#Trinidad_Head
#Tudor_Hill
#Tutuila

#location_list = ['Arrival_Heights','Barrow','Lauder','Mace_Head','Mauna_Loa','Niwot_Ridge','Ragged_Point','South_Pole','Trinidad_Head','Tudor_Hill','Tutuila']
location_list=['Mace_Head']
#location_list=['Cape_Verde']

#Model Version
#v90102
#v90103
mversion = 'v90103'

#Resolution
#4x5
#2x2.5
#0.5x0.666 Nested Europe
res = '4x5'

#Met.
#GEOS 5
#MERRA
met = 'GEOS 5'

#test = autojit(fast_lomb_phase.lomb)

for location in location_list:

	for species in species_list:
		units, obs_data_name, unit_cut, species_type, actual_species_name, obs_switch, model_cut_switch, ofac = modules.obs_variable_finder(species)

#Get variables for location chosen
		obsfile,loc_label,model_index = modules.location_check(location)

#set GAW_switch on or off. 'y' = multiple location GAW sim, 'n' = 1 location output
		GAW_switch = 'n'

# Read in the model output
		if GAW_switch == 'n':
			model , names = modules.readfile_GAW("binary_logs/GEOS_v90103_2x2.5_GAW_O3_logs.npy",model_index) #model index represents gaw location
		else:
			model , names = modules.readfile("binary_logs/GEOS_v90103_4x5_CV_logs.npy","001") #001 represents single location

# Processes the model date 
		date = model[:,0]
		time = model[:,1]
		model_time = modules.date_process(date,time)

#Define sampling intervals
		samp_spacing = 1./24.

#Convert model time array into numpy array
		model_time=np.array(model_time)

		counter = 0

# Remove posiibility of trying to plot species other than O3 at sites not Cape Verde
		if location != 'Cape_Verde':
			species_list = ['O3']

#set plotting area & background to white
		fig=plt.figure(figsize=(20,12))
		fig.patch.set_facecolor('white')
		ax = fig.add_subplot(1,1,1)

#Process Observations and date and normalise
		if location == 'Mace_Head':
			date, time, vals = modules.NOAA_data_reader_mace_head(obsfile)
			valid = vals >= 0
			#window = signal.hann(len(vals), sym=False)
			obs_vals = vals[valid]
			#window = window[valid]
			date = date[valid]
			time = time[valid]
			obs_time = modules.date_process(date,time)
		elif location == 'Cape_Verde':
			time, vals =  modules.read_CV(obsfile,obs_data_name,obs_switch)	
			valid = vals >= 0
			obs_vals = vals[valid]
			obs_time = time[valid]
		else:
			date, time, vals = modules.NOAA_data_reader(glob.glob(obsfile))
			valid = vals >= 0
			obs_vals = vals[valid]
			date = date[valid]
			time = time[valid]
			obs_time = modules.date_process(date,time)

		print type(vals)

	#standard_deviation_obs_p = np.std(obs_vals)
	#mean_obs_p = np.mean(obs_vals)
	#normal_obs = obs_vals-mean_obs_p
	#normal_obs = normal_obs/standard_deviation_obs_p

#Calculate variance of pre-processed obs data- should be 1 if normal
#standard_dev_obs = np.std(normal_var_2, dtype=np.float64)
#variance_obs = standard_dev_obs**2
#print 'Variance - pre-processed obs data= ', variance_obs

#Convert obs time array into numpy array
		obs_time=np.array(obs_time)

#Need to cut and normalise model data 
		if model_cut_switch == 0:
			k=names.index(species)
			model_cut = model[:,k]*unit_cut
		if model_cut_switch == 1:
			model_cut = modules.read_diff_species_wind(names,model)

	#standard_deviation_model_p = np.std(model_cut)
	#mean_model_p = np.mean(model_cut)
	#normal_model = model_cut-mean_model_p
	#normal_model = normal_model/standard_deviation_model_p

#Calculate variance of pre-processed model data- should be 1 if normal
#standard_dev_model = np.std(normal_model, dtype=np.float64)
#variance_model = standard_dev_model**2
#print 'Variance - pre-processed model data= ', variance_model

		np.savetxt('mho_obs_2006_2012',obs_vals)
		np.savetxt('mho_gc2x2.5_2006_2012',model_cut)

#align periods
		obs_time,model_time,obs_vals,model_cut = modules.align_periods(obs_time,model_time,obs_vals,model_cut)

		
		
#Define sampling frequency
		samp_freq = 24

		pi2 = np.pi*2
		daily_amplitude = 5
		hourly_step = 1./24
		syn_time = np.arange(0,365.25*3,hourly_step)
		waveform1 = daily_amplitude*(np.cos((pi2*syn_time/1)-(0)))
		waveform2 = daily_amplitude*(np.cos((pi2*syn_time/(365.25/2))-(0)))
		waveform3 = daily_amplitude*(np.cos((pi2*syn_time/365.25)-(0)))

		waveform = 50 + (waveform1+waveform2+waveform3)

		#test = (obs_time > 365) & (obs_time < 730) 
		#obs_time = obs_time[test]
		#obs_vals = obs_vals[test]


#Obs. lomb
		window = signal.hamming(len(obs_vals))
		#window = np.blackman(len(obs_vals))
		#window = np.kaiser(len(obs_vals),4)
		#window = signal.flattop(len(obs_vals))
		obs_mean = np.mean(obs_vals)
		obs_vals = obs_vals - obs_mean
		obs_vals = obs_vals*window
		#fa, fb, obs_amp, ph= lomb.fasper(obs_time,obs_vals)
		#annual_index = np.argmin(np.abs(1./fa - 365.25))
		#print 1./(fa[annual_index])
		#obs_periods = 1./fa
		NOUT = 0.5*4*1*len(obs_vals)
		NOUT = int(NOUT)
        #take lomb
		print NOUT
		print len(obs_time)
		print len(obs_vals)
		fa, fb, mag, ph, sig= lomb_phase.lomb(obs_time,obs_vals,NOUT)
		print sig
		amp_corr = 1./(sum(window)/len(window))
		obs_amp = mag * amp_corr

		#annual_index = np.argmin(np.abs(obs_periods - 1))
		#print 'annual obs= ', obs_amp[annual_index]
	#sigs = lomb.getSignificance(fa,fb,nout,ofac)    
#Divide output by sampling frequency
	#fb = fb/samp_freq
	#print len(fa)
	#print len(sigs)

	#signi_99, signi_95, signi_90, signi_50 = lomb.false_alarm_threshold_first_approx(nout)	
    
	#print signi_95
	#plt.axhline(y=signi_99, color = 'green', linestyle = '--')
	#plt.axhline(y=signi_95, color = 'green', linestyle = '--')
	#plt.axhline(y=signi_90, color = 'green', linestyle = '--')
	#plt.axhline(y=signi_50, color = 'green', linestyle = '--')
	#xy = (0.02,200)
	#plt.annotate('Obs. 95% Significance', xy, xytext=None, color='green')

		

#plot up
		#plt.loglog(1./fa, obs_amp, color = 'black', marker = 'x',markersize=2, label='Observational %s ' % (actual_species_name))
		#plt.plot(obs_time, obs_vals)
		#plt.loglog(1./fa_ph, amp_ph, color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='Obs. %s Lomb-Scargle' % (actual_species_name))
#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
#nyquist_freq_lomb_model = frequencies[-1]
#Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
#print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 



#Model 
	#fft
		#fft_array = fft(model_cut)
		#fft_mag = np.abs(fft_array)
		#fft_phase = [phase(i) for i in fft_array]
		#fft_phase = np.array(fft_phase)

		#fft_freqs = fftfreq(len(model_cut),1./24) 
		#valid = fft_freqs > 0
		#fft_mag, fft_phase, fft_freqs = fft_mag[valid], fft_phase[valid], fft_freqs[valid]

		#fft_mag = fft_mag/(len(obs_time)/2)
		#fft_periods = 1./fft_freqs

		#print len(fft_periods)
		#print fft_periods
		
		#plt.semilogx(fft_periods, fft_mag, color = 'red', marker = 'x', alpha = 0.75,markersize=2, label='GEOS-Chem %s %s %s %s FFT ' %(mversion,res,met,actual_species_name))
		#plt.plot(syn_time,waveform)

		daily_ave = []

		start = 0
		end = 31
		m_vals = [31,28,31,30,31,30,31,31,30,31,30,31,0]

		for i in range(12):
			daily_ave.append(np.mean(model_cut[start*24:end*24:8760]))
			start+=m_vals[i]
			end+=m_vals[i+1]

		
		#daily_ave_e = daily_ave[19:24] + daily_ave[0:19]
		plt.plot(range(12),daily_ave)
		#annual_index = np.argmin(np.abs(model_periods - 1))
		#print 'annual model phase= ', fft_phase[annual_index]
		

#Divide output by sampling freq
		window = signal.flattop(len(model_cut))
		#window = np.kaiser(len(model_cut),4)	
		model_mean = np.mean(model_cut)
		model_cut = model_cut - model_mean
		model_cut = model_cut*window
		#model_cut = model_cut[:-2]
		#model_time = model_time[:-2]
		#fc, fd, model_amp, model_ph= lomb.fasper(model_time,model_cut)
	
		#fft_array = fft(model_cut)
		#fft_mag = np.abs(fft_array)
		#fft_phase = [phase(i) for i in fft_array]
		#fft_phase = np.array(fft_phase)
		#fft_freqs = fftfreq(len(model_cut),hourly_step)

		#valid = fft_freqs > 0
		#fft_mag, fft_phase, fft_freqs = fft_mag[valid], fft_phase[valid], fft_freqs[valid]

		#model_periods = 1./fft_freqs

		NOUT = 0.5*4*1*len(model_cut)
		NOUT = int(NOUT)

		#fa, fb, fft_model, fft_phase= lomb_phase.lomb(model_time,model_cut,NOUT)
		#model_periods = 1./fa
		#amp_corr = 1./(sum(window)/len(window))
		#fft_model = fft_model * amp_corr	
		

		#model_periods = 1./fc
		#amp_corr = 1./(sum(window)/len(window))
		#model_amp = model_amp * amp_corr
		
	#fy = fy/samp_freq

	#signi_99, signi_95, signi_90, signi_50 = lomb.false_alarm_threshold_first_approx(nout)
	
	#print signi_95
	#plt.axhline(y=signi_99, color = 'blue', linestyle = '--')
	#plt.axhline(y=signi_95, color = 'blue', linestyle = '--')
	#plt.axhline(y=signi_90, color = 'blue', linestyle = '--')
	#plt.axhline(y=signi_50, color = 'blue', linestyle = '--')
	#xy = (0.02,100)
	#plt.annotate('Model 95% Significance', xy, xytext=None, color='blue')

#plot up
		#plt.plot(1./fc, model_amp, color = 'red', marker = 'x', alpha = 0.75,markersize=2, label='GEOS-Chem %s %s %s %s ' %(mversion,res,met,actual_species_name))

#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
#nyquist_freq_lomb_model = frequencies[-1]
#Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
#print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 
		def form2(x, pos):
			""" This function returns a string with 3 decimal places, given the input x"""
			return '%.2f' % x	

		def form5(x, pos):
			""" This function returns a string with 3 decimal places, given the input x"""
			return '%.6f' % x

		xformatter = FuncFormatter(form2)
		yformatter = FuncFormatter(form5)
		#plt.xlim(1e-2,10000)
		#plt.ylim(1e-6,1000)
		plt.xlim(0,11)
		#plt.ylim(0.1,1.3)
		plt.grid(True)
		leg=plt.legend(loc=4, prop={'size':21})
		#leg.get_frame().set_alpha(0.4)
		#ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
		#ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
		#plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
		#plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))
		plt.xlabel('Day',fontsize = 21)
		#plt.xlabel('Period (Days)', fontsize=21)
		#plt.ylabel('Magnitude (%s)'%units, fontsize=21)
		#plt.title('Spectral Magntitudes - %s'%loc_label,fontsize=21)
		#plt.tick_params(axis='y',which='both',bottom='off',top='off',labelbottom='off')
		plt.xticks([0,1,2,3,4,5,6,7,8,9,10,11],['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'])
		plt.gca().yaxis.set_major_locator(plt.NullLocator())
		plt.title('Annual Continental O3', fontsize = 21)
		counter+=1


		#plt.savefig('plots/v90103_4x5_MERRA/normal/%s_%s_%s_%s.png'%(actual_species_name,location,mversion,res),dpi=fig.dpi)

		plt.show()










