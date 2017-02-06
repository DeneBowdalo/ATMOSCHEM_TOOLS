import numpy as np
import csv
import glob
import datetime as datetime
import matplotlib.pyplot as plt
import logging as log
#import nappy
from cmath import *
import scipy.stats as stats
import lomb_phase
import modules
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq, rfft
import matplotlib.ticker
from matplotlib.ticker import FuncFormatter
from scipy import signal
from rpy2.robjects import IntVector, Formula
from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri as rpyn
import rpy2.robjects as robjects


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
location_list = ['Mace_Head']

#Model Version
#v90102
#v90103
mversion = 'v90103'

#Resolution
#4x5
#2x2.5
#0.5x0.666 Nested Europe
res = '2x2.5'

#Met.
#GEOS 5
#MERRA
met = 'GEOS 5'

for location in location_list:

# Remove posiibility of truing to plot species other than O3 at sites not Cape Verde
	if location != 'Cape_Verde':
		species_list = ['O3']

#Get variables for location chosen
	obsfile,loc_label,model_index = modules.location_check(location)

	for species in species_list:
		units, obs_data_name, unit_cut, species_type, actual_species_name, obs_switch, model_cut_switch, ofac = modules.obs_variable_finder(species)

#set GAW_switch on or off. 'y' = multiple location GAW sim, 'n' = 1 location output
		GAW_switch = 'y'

# Read in the model output
		if GAW_switch == 'y':
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

		for species in species_list:
			units, obs_data_name, unit_cut, species_type, actual_species_name, obs_switch, model_cut_switch, ofac = modules.obs_variable_finder(species)

#set plotting area & background to white
		fig=plt.figure(figsize=(20,12))
		fig.patch.set_facecolor('white')
		ax = fig.add_subplot(1,1,1)

#Process Observations and date and normalise
		if location == 'Mace_Head':
			date, time, vals = modules.NOAA_data_reader_mace_head(obsfile)
			valid = vals >= 0
			obs_vals = vals[valid]
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

#align periods
		obs_time,model_time,obs_vals,model_cut = modules.align_periods(obs_time,model_time,obs_vals,model_cut)

#Define sampling frequency
		samp_freq = 24

#Obs. lomb
		window = signal.hamming(len(obs_vals))
		obs_mean = np.mean(obs_vals)
		obs_vals = obs_vals - obs_mean
		obs_vals = obs_vals*window

		#processing to determine freq range
		OFAC = 1
		HIFAC = 1
								
		NOUT = 0.5*OFAC*HIFAC*len(obs_vals)
		NOUT = int(NOUT)
		SAMP_R = 1./24
		
		#take lomb
		fa, fb, mag, ph, sig, i_freq, s_val_ind, end_val_ind= lomb_phase.lomb(obs_time,obs_vals,NOUT,OFAC,SAMP_R)
		#limit arrays by valid indices determine in fortran code
		fa = fa[s_val_ind:end_val_ind]
		fb = fb[s_val_ind:end_val_ind]
		mag = mag[s_val_ind:end_val_ind]
		ph = ph[s_val_ind:end_val_ind]
		sig = sig[s_val_ind:end_val_ind]
																		
		obs_periods = 1./fa

		obs_periods = np.log(obs_periods)
		mag = np.log(mag)
		
		reversed_obs_periods = obs_periods[::-1]
		reversed_obs_amp = mag[::-1]
		reversed_sig = sig[::-1]

		a = sorted(sig)
		print a[0:300]

		sig_test = sig >= 99.9
		#reversed_sig = reversed_sig[sig_test] 
		#obs_periods_c = obs_periods[sig_test]
		#fb_c = fb[sig_test]

		un_sig_indices = np.nonzero(sig_test == False)
		sig_indices = np.nonzero(sig_test == True)

		print obs_periods[sig_indices]

		un_sig_indices = np.array(un_sig_indices)
		un_sig_indices = un_sig_indices.flatten()
		
		fa_c = fa[un_sig_indices]
		fb_c = fb[un_sig_indices]
		mag_c = mag[un_sig_indices]
		ph_c = ph[un_sig_indices]
		sig_c = sig[un_sig_indices]
		obs_periods_c = obs_periods[un_sig_indices]

		reversed_obs_periods = obs_periods_c[::-1]
		reversed_obs_amp = mag_c[::-1]
		normal_obs_periods = obs_periods_c
		normal_obs_amp = mag_c
		#reversed_sig = sig[::-1]

		#print len(reversed_obs_amp)

		print len(obs_periods)
		print len(obs_periods_c)

		#sig_periods, sig_mag, sig_indices = modules.peak_rank(reversed_obs_periods,reversed_obs_amp)
		#print type(sig_indices)
#Get significant periods
		#reversed_obs_periods = obs_periods[::-1]
		#reversed_obs_amp = obs_amp[::-1]
		#sig_periods, sig_mag, sig_indices = modules.significant_peaks(reversed_obs_periods,reversed_obs_amp)

		#obs_amp = np.log(obs_amp)

                #limit array to between 0 and 50.
		test = normal_obs_periods < 50
		normal_obs_periods = normal_obs_periods[test]
		normal_obs_amp = normal_obs_amp[test]
		#test = normal_obs_periods > 1
		#normal_obs_periods = normal_obs_periods[test]
		#normal_obs_amp = normal_obs_amp[test]

#Smooth using exponential smoother, last 2 numbers inputted are: smoothing window & cut-off point
		#cut_obs_periods, smoothed_obs = modules.ema_keep_peaks(reversed_obs_amp,reversed_obs_periods,40,1200000)
		#cut_obs_periods, smoothed_obs = modules.just_smooth(reversed_obs_amp,reversed_obs_periods,30000)
		cut_obs_periods, smoothed_obs = modules.log_bin_smooth(len(normal_obs_periods),20,normal_obs_periods,normal_obs_amp)
		#cut_obs_periods, smoothed_obs = modules.automatic_smooth_2(reversed_obs_amp,reversed_obs_periods,100,sig_indices)
		#smoothed_obs = np.exp(smoothed_obs)

		#limit array to between 0 and 50.
		#test = cut_obs_periods < 30
		#cut_obs_periods = cut_obs_periods[test]
		#smoothed_obs = smoothed_obs[test]
		#test = cut_obs_periods > 1
		#cut_obs_periods = cut_obs_periods[test]
		#smoothed_obs = smoothed_obs[test]

#plot up
		#plt.loglog(obs_periods, mag, color = 'black', marker = 'x',markersize=2, label = 'normal')
		#plt.loglog(obs_periods_c, mag_c, color = 'red', marker = 'x', alpha = 0.2,markersize=2, label = 'cut')
		plt.loglog(cut_obs_periods, smoothed_obs, color = 'black', marker = 'x', alpha = 1,markersize=2, label='Observational %s ' % (actual_species_name))
        #plt.plot(obs_time,obs_vals, color = 'black', marker = 'x', markersize = 2, label='Observational %s ' % (actual_species_name))
		
		#find breakpoint between turbulence and synoptic systems
		cut_obs_periods = np.array(cut_obs_periods)
		smoothed_obs = np.array(smoothed_obs)
		seg_mod = importr("segmented")
		stats = importr('stats')
		r = robjects.r
		#convert to r numerical vectors
		cut_obs_periods=rpyn.numpy2ri(cut_obs_periods)
		smoothed_obs=rpyn.numpy2ri(smoothed_obs)
		robjects.globalenv["x"] = cut_obs_periods
		robjects.globalenv["y"] = smoothed_obs
		fit = r.lm("y ~ x")

		#fit linear model to data
		x_fmla = Formula('~ x')

		#use segmented module
		#seg_out = seg_mod.segmented(fit, seg_Z = x_fmla, psi = robjects.IntVector((1,10)))
		seg_out = seg_mod.segmented(fit, seg_Z = x_fmla, psi = 10)
		print seg_mod.slope(seg_out)
		breakpoints = seg_out.rx('psi')
		breakpoints = rpyn.ri2numpy(breakpoints)
		breakpoint_1 = breakpoints[0,1].tolist()
		#breakpoint_2 = breakpoints[1,1].tolist()
		print 'breakpoint = ',  breakpoints[0,1]#, breakpoints[1,1]
		#plt.axvline(x=1, color = 'red', linestyle = '--')
		plt.axvline(x=breakpoint_1, color = 'blue', linestyle = '--')
		#plt.axvline(x=breakpoint_2, color = 'blue', linestyle = '--')
#Model lomb
		#model_cut = model_cut - np.mean(model_cut)
		#fc, fd, model_amp,ph = lomb.fasper(model_time,model_cut)
		#model_periods = 1./fc
		
#Reverse array for smoothing
		#reversed_model_amp = model_amp[::-1]
		#reversed_model_periods = model_periods[::-1]
		#sig_periods, sig_mag, sig_indices = modules.peak_rank(reversed_model_periods,reversed_model_amp)	
	#sig_periods, sig_mag, sig_indices = modules.significant_peaks(reversed_model_periods,reversed_model_amp)
		#model_amp = np.log(model_amp)
		#reversed_model_periods = model_periods[::-1]
		#reversed_model_amp = model_amp[::-1]

#Smooth using exponetial smoother, last 2 numbers inputted are: smoothing window & cut-off point
		#cut_model_periods, smoothed_model = modules.ema_keep_peaks(reversed_model_amp,reversed_model_periods,40,1200000)
		#cut_model_periods, smoothed_model = modules.automatic_smooth_2(reversed_model_amp,reversed_model_periods,100,sig_indices)
		#smoothed_model = np.exp(smoothed_model)

#plot up
		#plt.loglog(model_periods,model_amp, color = 'red', marker = 'x', alpha = 0.2,markersize=2)
		#plt.plot(cut_model_periods, smoothed_model, color = 'red', marker = 'x', alpha = 1,markersize=2, label='GEOS-Chem %s %s %s %s '%(mversion,res,met,actual_species_name))
        #plt.plot(model_time,model_cut, color = 'red', marker = 'x', markersize = 2, label='GEOS-Chem %s %s %s %s '%(mversion,res,met,actual_species_name)) 

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

		#fft_mag = fft_mag/(len(model_time)/2)

#Divide output by sampling frequency
		#fft_mag = np.log(fft_mag)

#Reverse array for smoothing
		#reversed_fft_mag = fft_mag[::-1]
		#fft_periods = 1./fft_freqs
		#reversed_fft_periods = fft_periods[::-1]
        #protect_periods = [1,345]

#Smooth using exponetial smoother, last 2 numbers inputted are: smoothing window & cut-off point
		#cut_fft_periods, smoothed_fft = modules.ema_keep_peaks(reversed_fft_mag,reversed_fft_periods,40,1200000)
		#smoothed_fft = np.exp(smoothed_fft)

		#plt.loglog(cut_fft_periods, smoothed_fft, color = 'red', marker = 'x', alpha = 0.75,markersize=2, label='GEOS-Chem %s %s %s %s FFT ' %(mversion,res,met,actual_species_name))

    #lomb

	#fx, fy, nout, jmax, prob2 = lomb.fasper(model_time,normal_model, ofac, 1.)
#Divide output by sampling frequency
	#fy = fy/samp_freq
	#fy = np.log(fy)

#Reverse array for smoothing
	#reversed_fy = fy[::-1]
	#model_periods = 1./fx
	#reversed_model_periods = model_periods[::-1]

#Smooth using exponetial smoother, last 2 numbers inputted are: smoothing window & cut-off point
	#cut_model_periods, smoothed_model = modules.ema(reversed_fy,reversed_model_periods,10,120)
	#smoothed_model = np.exp(smoothed_model)

#plot up
	#plt.loglog(cut_model_periods, smoothed_model, color = 'red', marker = 'x', alpha = 0.75,markersize=2, label='GEOS-Chem %s %s %s %s Smoothed ' %(mversion,res,met,actual_species_name))

#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
#nyquist_freq_lomb_model = frequencies[-1]
#Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
#print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 

		def form2(x, pos):
			""" This function returns a string with 3 decimal places, given the input x"""
			return '%.3f' % x

		def form5(x, pos):
			""" This function returns a string with 3 decimal places, given the input x"""
			return '%.5f' % x

		xformatter = FuncFormatter(form2)
		yformatter = FuncFormatter(form2)
		#plt.xlim(0.998,1.002)
        #plt.xticks((0,365,730,1095,1460,1825,2190),('2006','2007','2008','2009','2010','2011','2012'))
		#plt.xticks((0.998,1.000,1.002))
		#plt.xlim(1e-2,10000)
		#plt.ylim(1e-5,10)
		#plt.xlim(0,20)
		plt.grid(True)
		leg=plt.legend(loc=2, prop={'size':21})
		leg.get_frame().set_alpha(0.4)
		ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
		ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
		plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
		plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))
		plt.xlabel('Period (Days)', fontsize=21)
		plt.ylabel('Magnitude (%s)'%units, fontsize=21)
        #plt.xlabel('Time (Year)', fontsize = 21)
        #plt.ylabel('Concentration (ppbv)', fontsize = 21)
		plt.title('Smoothed Spectral Magntitudes determined by Lomb-Scargle Periodogram - %s'%loc_label,fontsize=21)
        #plt.title('Time Series at %s'%loc_label,fontsize = 21)
		counter+=1
	
		#plt.savefig('plots/v90103_4x5_MERRA/smoothed/%s_%s_%s_%s.png'%(actual_species_name,location,mversion,res),dpi=fig.dpi)

		plt.show()










