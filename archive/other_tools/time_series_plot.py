import numpy as np
import csv
import glob
import datetime as datetime
import matplotlib.pyplot as plt
import logging as log
#import nappy
import scipy.stats as stats
import lomb
import modules
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
obsfile,loc_label,model_index = modules.location_check(location)

#set GAW_switch on or off. 'y' = multiple location GAW sim, 'n' = 1 location output
GAW_switch = 'n'

# Read in the model output
if GAW_switch == 'y':
	model , names = modules.readfile_GAW("binary_logs/GEOS_v90103_4x5_GAW_O3_logs.npy",model_index) #model index represents gaw location
else:
	model , names = modules.readfile("binary_logs/GEOS_v90103_4x5_CV_logs.npy","001") #001 represents single location

# Processes the model date 
date = model[:,0]
time = model[:,1]
model_time = modules.date_process(date,time)

print model_time 

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
		valid = vals > 0
		obs_vals = vals[valid]
		obs_time = time[valid]
		interp_vals = np.interp(time,obs_time,obs_vals)	

	else:
		date, time, vals = modules.NOAA_data_reader(glob.glob(obsfile))
		valid = vals >= 0
		obs_vals = vals[valid]
		date = date[valid]
		time = time[valid]
		obs_time = modules.date_process(date,time)

#Convert obs time array into numpy array
	obs_time=np.array(obs_time)

#Need to cut and normalise model data 
	if model_cut_switch == 0:
		k=names.index(species)
		model_cut = model[:,k]*unit_cut
	if model_cut_switch == 1:
		model_cut = modules.read_diff_species_wind(names,model)

#cut obs to year 2007
print len(time)
print len(interp_vals)


test_2007 = (time >= 730) & (time < 1096)
obs_vals = interp_vals[test_2007]
obs_time  = time[test_2007]

obs_time = range(len(obs_time))


#cut model to year 2010
test_2010 = (model_time >= 1461) & (model_time < 1826)
model_cut = model_cut[test_2010]
time_2010 = range(len(model_cut))

year_time = np.arange(0,365,1./24)
print len(year_time)
print len(obs_vals)

#plot up
plt.plot(obs_time, obs_vals, color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='2008 %s Obs. %s' % (loc_label,actual_species_name))
plt.plot(time_2010, model_cut, color = 'red', marker = 'x', alpha = 0.75,markersize=2, label='2010 GEOS-Chem %s %s %s %s ' %(mversion,res,met,actual_species_name))

#plt.xlim(1e-2,1e4)
plt.ylim(0,60)
plt.tick_params(axis='x', which='major', labelsize=14)
plt.tick_params(axis='y', which='major', labelsize=14)
plt.tick_params(direction='out', pad=0.00005)
plt.grid(True)
leg=plt.legend(loc=4, prop={'size':21})
leg.get_frame().set_alpha(0.4)
plt.xlabel('Period (Hours)', fontsize=16, labelpad = 8)
plt.ylabel('%s (%s)'%(species_type,units), fontsize=16, labelpad = 8)
plt.title('Time Series',fontsize=22)
counter+=1

	#plt.savefig('plots/v90103_nestedeurope_1year_GEOS5/normal/%s_%s_%s_%s.png'%(actual_species_name,location,mversion,res),dpi=fig.dpi)

plt.show()










