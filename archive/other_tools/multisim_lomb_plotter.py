import numpy as np
import csv
import glob
import datetime as datetime
import matplotlib.pyplot as plt
from matplotlib import cm
import logging as log
import nappy
import scipy.stats as stats
import lomb
import modules
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq, rfft
#select speci (or species) & location. Only Cape Verde has multiple species

#O3
#CO
#NO
#NO2
#C2H6
#C3H8
#DMS
#TRA_6
#ACET
#GMAO_TEMP
#GMAO_PSFC
#GMAO_WIND
#GMAO_RADSW
#GMAO_RHUM

species = 'O3' 

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

location = 'Mace_Head'

# Remove posiibility of truing to plot species other than O3 at sites not Cape Verde
if location != 'Cape_Verde':
    species_list = ['O3']

#Model Version
#v90102
#v90103
mversion = ['v90103','v90103','v90103','v90103']

#Resolution
#4x5
#2x2.5
#0.5x0.666 Nested Europe
res = ['4x5','4x5','2x2.5','0.5x0.666 Nested Europe']

#Met.
#GEOS 5
#MERRA
met = ['GEOS 5','MERRA','GEOS 5','GEOS 5']

#Get variables for location chosen
obsfile,loc_label,model_index = modules.location_check(location)

#set GAW_switch on or off. 'y' = multiple location GAW sim, 'n' = 1 location output
GAW_switch = ['y','n','y','y']

counter = 0

#set plotting area & background to white
fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')

#species details
units, obs_data_name, unit_cut, species_type, actual_species_name, obs_switch, model_cut_switch, ofac = modules.obs_variable_finder(species)

#Process Observations and date and normalise
if location == 'Mace_Head':
	date, time, vals = modules.NOAA_data_reader_mace_head(obsfile)
	valid = vals >= 0
	obs_vals = vals[valid]
	date = date[valid]
	time = time[valid]
	obs_time = modules.date_process(date,time)
elif location == 'Cape_Verde':
	time, vals =  modules.CV_data_reader(obsfile,obs_data_name,obs_switch)
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

standard_deviation_obs_p = np.std(obs_vals)
mean_obs_p = np.mean(obs_vals)
normal_obs = obs_vals-mean_obs_p
normal_obs = normal_obs/standard_deviation_obs_p

#Calculate variance of pre-processed obs data- should be 1 if normal
#standard_dev_obs = np.std(normal_var_2, dtype=np.float64)
#variance_obs = standard_dev_obs**2
#print 'Variance - pre-processed obs data= ', variance_obs

#Convert obs time array into numpy array
obs_time=np.array(obs_time)

#Define sampling frequency
samp_freq = 24

#Obs. lomb
fa, fb, nout, jmax, prob2 = lomb.fasper(obs_time,normal_obs, ofac, samp_freq)
#Divide output by sampling frequency
fb = fb/samp_freq
fb = np.log(fb)

#Reverse array for smoothing
reversed_fb = fb[::-1]
obs_periods = 1./fa
reversed_obs_periods = obs_periods[::-1]

#Smooth using exponetial smoother, last 2 numbers inputted are: smoothing window & cut-off point
cut_obs_periods, smoothed_obs = modules.ema_keep_peaks(reversed_fb,reversed_obs_periods,40,12000000)
smoothed_obs = np.exp(smoothed_obs)

#plot up
plt.loglog(cut_obs_periods, smoothed_obs, color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='%s Obs. %s Smoothed ' % (loc_label,actual_species_name))

#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
#nyquist_freq_lomb_model = frequencies[-1]
#Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
#print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 


#loop through diff models and get respective obs. data
for i in range(len(mversion)):

	if GAW_switch[counter] == 'n':
		if mversion[counter] == 'v90102':
			if res[counter] == '4x5':
				if met[counter] == 'GEOS 5':
					model , names = modules.readfile("binary_logs/GEOS_v90102_4x5_CV_logs.npy","001") #001 represents single location
				if met[counter] == 'MERRA':
					model , names = modules.readfile("","001") #001 represents single location
			if res[counter] == '2x2.5':
				if met[counter] == 'GEOS 5':
					model , names = modules.readfile("binary_logs/GEOS_v90103_2x25_CV_logs.npy","001") #001 represents single location
				if met[counter] == 'MERRA':
					model , names = modules.readfile("","001") #001 represents single location

 
		if mversion[counter] == 'v90103':
			if res[counter] == '4x5':
				if met[counter] == 'GEOS 5':
					model , names = modules.readfile("binary_logs/GEOS_v90103_4x5_CV_logs.npy","001") #001 represents single location
				if met[counter] == 'MERRA':
					model , names = modules.readfile("binary_logs/GEOS_v90103_4x5_MERRA_Mace_Head_logs.npy","001") #001 represents single location
			if res[counter] == '2x2.5':
				if met[counter] == 'GEOS 5':
					model , names = modules.readfile("binary_logs/GEOS_v90103_2x25_CV_logs.npy","001") #001 represents single location
				if met[counter] == 'MERRA':
					model , names = modules.readfile("","001") #001 represents single location
			if res[counter] == '0.5x0.666 Nested Europe':
				if met[counter] == 'GEOS 5':
					model , names = modules.readfile("","001") #001 represents single location
				if met[counter] == 'MERRA':
					model , names = modules.readfile("","001") #001 represents single location


	if GAW_switch[counter] == 'y':
		if mversion[counter] == 'v90102':
			if res[counter] == '4x5':
				if met[counter] == 'GEOS 5':
					model , names = modules.readfile_GAW("",model_index) #model index represents gaw location                                     
				if met[counter] == 'MERRA':
					model , names = modules.readfile_GAW("",model_index) #model index represents gaw location 
			if res[counter] == '2x2.5':        
				if met[counter] == 'GEOS 5':
					model , names = modules.readfile_GAW("",model_index) #model index represents gaw location 
				if met[counter] == 'MERRA':
					model , names = modules.readfile_GAW("",model_index) #model index represents gaw location
  
		if mversion[counter] == 'v90103':
			if res[counter] == '4x5':
				if met[counter] == 'GEOS 5':
					model , names = modules.readfile_GAW("binary_logs/GEOS_v90103_4x5_GAW_O3_logs.npy",model_index) #model index represents gaw location                                     
				if met[counter] == 'MERRA':
					model , names = modules.readfile_GAW("",model_index) #model index represents gaw location
			if res[counter] == '2x2.5':        
				if met[counter] == 'GEOS 5':
					model , names = modules.readfile_GAW("binary_logs/GEOS_v90103_2x2.5_GAW_O3_logs.npy",model_index) #model index represents gaw location                                 
				if met[counter] == 'MERRA':
					model , names = modules.readfile_GAW("",model_index) #model index represents gaw location
			if res[counter] == '0.5x0.666 Nested Europe':
				if met[counter] == 'GEOS 5':
					model , names = modules.readfile("binary_logs/GEOS_v90103_nested_europe_GAW_logs.npy",model_index) #model index represents gaw location                                     
				if met[counter] == 'MERRA':
					model , names = modules.readfile_GAW("",model_index) #model index represents gaw location  
              
# Processes the model date 
	date = model[:,0]
	time = model[:,1]


		
	model_time = modules.date_process(date,time)

	print counter
	print res[counter]
	print met[counter]

	#print date	
	#print time
#Define sampling intervals
	samp_spacing = 1./24.

#Convert model time array into numpy array
	model_time=np.array(model_time)

#Need to cut and normalise model data 
	if model_cut_switch == 0:
		k=names.index(species)
		model_cut = model[:,k]*unit_cut
	if model_cut_switch == 1:
		model_cut = module.read_diff_species_wind()

	standard_deviation_model_p = np.std(model_cut)
	mean_model_p = np.mean(model_cut)
	normal_model = model_cut-mean_model_p
	normal_model = normal_model/standard_deviation_model_p

#Calculate variance of pre-processed model data- should be 1 if normal
#standard_dev_model = np.std(normal_model, dtype=np.float64)
#variance_model = standard_dev_model**2
#print 'Variance - pre-processed model data= ', variance_model

#Model lomb
	fx, fy, nout, jmax, prob2 = lomb.fasper(model_time,normal_model, ofac, samp_freq)
#Divide output by sampling frequency
	fy = fy/samp_freq
	fy = np.log(fy)

#Reverse array for smoothing
	reversed_fy = fy[::-1]
	model_periods = 1./fx
	reversed_model_periods = model_periods[::-1]

#Smooth using exponetial smoother, last 2 numbers inputted are: smoothing window & cut-off point
	cut_model_periods, smoothed_model = modules.ema_keep_peaks(reversed_fy,reversed_model_periods,40,1200000)
	smoothed_model = np.exp(smoothed_model)

#plot up
	plt.loglog(cut_model_periods, smoothed_model, color=cm.jet(1.*i/len(mversion)) , marker = 'x', alpha = 0.75,markersize=2, label='GEOS-Chem %s %s %s %s Smoothed ' %(mversion[counter],res[counter],met[counter],actual_species_name))

#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
#nyquist_freq_lomb_model = frequencies[-1]
#Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
#print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 

	plt.grid(True)
	leg=plt.legend(loc=4, prop={'size':16})
	leg.get_frame().set_alpha(0.4)
	plt.xlabel('Period (Days)', fontsize=21)
	plt.ylabel(r'PSD $(ppb^{2}/days^{-1})$', fontsize=21)
	plt.title('Smoothed Lomb-Scargle Periodogram.',fontsize=21)
	counter+=1
	#plt.xlim(0.07,10)
	#plt.ylim(1e-7,10)
plt.show()













