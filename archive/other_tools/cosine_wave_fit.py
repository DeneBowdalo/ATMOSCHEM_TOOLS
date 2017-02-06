import numpy as np
import csv
import glob
import datetime as datetime
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import logging as log
import nappy
import scipy.stats as stats
import lomb
import modules
import pylab
from pylab import *
from matplotlib.backends.backend_pdf import PdfPages
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq, rfft
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

location = 'Cape_Verde'

# Remove posiibility of truing to plot species other than O3 at sites not Cape Verde
if location != 'Cape_Verde':
	species_list = ['O3']

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

#Get variables for location chosen
#obsfile,loc_label,model_index = modules.location_check(location)
#obsfile,loc_label = modules.long_obs_data(location)

def date_process_long(date,time):
    year=(date//10000)
    month=((date-year*10000)//100)
    day=(date-year*10000-month*100)

    hour=time//100
    min=(time-hour*100)

    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(1973,1,1,0,0,0) \
              for i in range(len(year))]

    processed_dates=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
    return processed_dates

#set GAW_switch on or off. 'y' = multiple location GAW sim, 'n' = 1 location output
GAW_switch = 'n'

if GAW_switch == 'y':
	model , names = modules.readfile_GAW("binary_logs/GEOS_v90103_nested_europe_GAW_logs.npy",model_index) #model index represents gaw location
else:
	model , names = modules.readfile("binary_logs/GEOS_v90103_4x5_CV_logs.npy","001") #001 represents single location

# Processes the model date 
date = model[:,0]
time = model[:,1]

print model
#Define sampling intervals
#samp_spacing = 1./24.

#Convert model time array into numpy array

counter = 0

for species in species_list:
	units, obs_data_name, unit_cut, species_type, actual_species_name, obs_switch, model_cut_switch, ofac = modules.obs_variable_finder(species)

#set plotting area & background to white
	fig=plt.figure(figsize=(20,12))
	fig.patch.set_facecolor('white')
	ax = fig.add_subplot(111)
#Process Observations and date and normalise
	#if location == 'Mace_Head':
	#	date, time, vals = modules.NOAA_data_reader_mace_head(obsfile)
	#	valid = vals >= 0
	#	obs_vals = vals[valid]
	#	date = date[valid]
	#	time = time[valid]
	#	obs_time = modules.date_process(date,time)
	#elif location == 'Cape_Verde':
	#	time, vals =  modules.CV_data_reader(obsfile,obs_data_name,obs_switch)	
	#	valid = vals >= 0
	#	obs_vals = vals[valid]
	#	obs_time = time[valid]
	#else:
	#	date, time, vals = modules.NOAA_data_reader(glob.glob(obsfile))
	
	years = [str(i)[0:4] for i in date]
	years = map(int, years)
	years = np.array(years)
	current_year = years[0]
	final_year = years[-1]
	years_to_iterate = (final_year - current_year) + 1 
		
#Need to cut and normalise model data 
	if model_cut_switch == 0:
		k=names.index(species)
		model_cut = model[:,k]*unit_cut
	if model_cut_switch == 1:
		model_cut = modules.read_diff_species_wind(names,model)

	model_time = modules.date_process(date,time)
	model_time = np.array(model_time)

	ave_species = np.average(model_cut)

	start = 0
	end = 8760
	leap_year_count = 1
	timeoffsets = []

	for i in range(len(model_cut)//8766): # for 365.25 days		
		yearly_cut = model_cut[start:end]
		time_cut  = model_time[start:end]
		start = end
		end += 8760
		leap_year_count+=1

		#limit data to within interquartile range
		#upperQuartile_year = stats.scoreatpercentile(yearly_cut,75)
		#lowerQuartile_year = stats.scoreatpercentile(yearly_cut,25)
		#yearly_valid_cut_test = (yearly_cut >= lowerQuartile_year) & (yearly_cut <= upperQuartile_year)  
		#yearly_valid_cut = yearly_cut[yearly_valid_cut_test]		
		#print len(yearly_valid_cut)

		#Split yearly data into chunks
		new_yearly_cut_chunks = []
		yearly_cut = np.array(yearly_cut)
		yearly_cut_chunks = np.split(yearly_cut, 4)
		for i in yearly_cut_chunks:
			new_yearly_cut_chunks.append(np.average(i))
		maximum = new_yearly_cut_chunks.index(np.max(new_yearly_cut_chunks))
		yearly_amplitude_test = np.max(yearly_cut_chunks[maximum]) == yearly_cut
		
		#yearly_amplitude_test = np.max(yearly_cut) == yearly_cut
		yearly_timeoffset = time_cut[yearly_amplitude_test] 		
		yearly_timeoffset = np.average(yearly_timeoffset)
		timeoffsets.append(yearly_timeoffset)
		if leap_year_count == 4:
			end +=24

	n=0
	adjusted = []
	for i in timeoffsets:
		adjusted_offset = i - (365.25*n)	
		n+=1
		adjusted.append(adjusted_offset)
	time_offset = np.average(adjusted)
	
# Calculate the interquartile range, equal to amplitude                                                                     
	sorted_data = np.sort(model_cut)                                                                                                   
	upperQuartile = stats.scoreatpercentile(sorted_data,75)                                                                      
	lowerQuartile = stats.scoreatpercentile(sorted_data,25)                                                                      
	IQR = upperQuartile - lowerQuartile  	
	amplitude = IQR			
		
	print 'Amplitude = ', amplitude
	print 'Time offset = ', time_offset
	print 'Vertical_Average_Shift = ', ave_species
	print 'Period = ', 365.25
	#print 'Energy Spectral Density = ', amplitude**2
	
	#fit sine wave
	o3_wave = ave_species + (amplitude*np.cos(2*np.pi*(model_time-time_offset)/365.25))

	#print 'Power Spectral Density = ', (amplitude**2)/len(o3_wave)
	difference  = model_cut / o3_wave
	#print 'Ratio of Actual Wave to Idealised Cosine Wave = ', np.average(difference) 
	
# Remove annual period
	

#Remove half- annual signal


#Remove diurnal signal

	samp_spacing=float(1./24.)
	fft_total=fft(o3_wave)
	print len(fft_total)
	print fft_total[1]
	fft_model=np.abs(fft(o3_wave))#**2 #calculate the magnitude of signal components
	frequencies = fftfreq(o3_wave.size, d = samp_spacing)
	keep_freq = frequencies>0
	frequencies, fft_model = frequencies[keep_freq], fft_model[keep_freq]
	model_periods = 1 /frequencies
	print fft_model[0]
	fft_model = fft_model/len(model_cut)
	
	#print len(fft_model)
	#print fft_model[1]
	#annual_amplitude_test =(364 <model_periods) & (model_periods < 366)
	#fft_annual_amplitude = fft_model[annual_amplitude_test]
	#print fft_annual_amplitude
	#print 'FFT Annual Amplitude = ', fft_annual_amplitude  
	#print 'FFT Annual Amplitude corrected = ', fft_annual_amplitude*2

	#ax1.loglog(model_periods, fft_model , color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='FFT %s' % (actual_species_name))
	

	#Plot up
	plt.plot(model_time,model_cut,marker ='x', markersize = 1, linestyle= 'None',color='red', label='GEOS-Chem %s %s %s %s ' %(mversion,res,met,actual_species_name)) 
	plt.plot(model_time,o3_wave,color='blue', label='Idealised Annual Cosine wave fit')	
	plt.grid(True)
	leg=plt.legend(loc=4, prop={'size':21})
	leg.get_frame().set_alpha(0.4)
	plt.xlabel('Period (Days)', fontsize=21)
	plt.ylabel('%s (%s)'%(actual_species_name,units), fontsize=21)
	plt.title('Time series Cosine wave fit',fontsize=21)
	#ax.text(0, 0.99,'Ratio of Actual Wave to Idealised Cosine Wave = %10.7f' %(np.average(difference)), ha='left', va='center', transform=ax.transAxes, fontweight ='bold')

	

#plt.ylim(-4,4)
#plt.xlim(0,365)
	plt.show()


print len(fft_model)
#print fft_model[1]
annual_amplitude_test =(364 <model_periods) & (model_periods < 366)
fft_annual_amplitude = fft_model[annual_amplitude_test]
print fft_annual_amplitude
print 'FFT Annual Amplitude = ', fft_annual_amplitude
print 'FFT Annual Amplitude corrected = ', fft_annual_amplitude*2

plt.loglog(model_periods, fft_model , color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='FFT %s' % (actual_species_name))

plt.show()

#set plotting area & background to white
#fig=plt.figure(figsize=(20,12))
#fig.patch.set_facecolor('white')

#samp_freq = 24 

#Model lomb
#fx, fy, nout, jmax, prob2 = lomb.fasper(model_time,o3_wave, ofac, samp_freq)
#Divide output by sampling frequency
#fy = fy/samp_freq

#signi_99, signi_95, signi_90, signi_50 = lomb.false_alarm_threshold_first_approx(nout)
#plt.axhline(y=signi_95, color = 'blue', linestyle = '--')


#plt.loglog(1./fx, fy, color = 'red', marker = 'x', alpha = 0.75,markersize=2, label='Fitted Cosine')
#plt.xlim(1e-2,1e4)
    #plt.ylim(1e-300,1e10)
#plt.grid(True)
#leg=plt.legend(loc=4, prop={'size':21})
#leg.get_frame().set_alpha(0.4)
#plt.xlabel('Period (Days)', fontsize=21)
#plt.ylabel(r'PSD $(ppb^{2}/days^{-1})$', fontsize=21)
#plt.title('Lomb-Scargle Periodogram',fontsize=21)
#plt.show()







