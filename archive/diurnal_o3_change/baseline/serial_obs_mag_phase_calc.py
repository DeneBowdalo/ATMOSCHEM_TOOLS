import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from cmath import *
from math import *
import csv
import datetime
import lomb_phase
from scipy import signal
import multiprocessing
import datetime
import time

start = datetime.datetime.now()

daily_period = 1
half_annual_period = (365.25/2)
annual_period = 365.25

#results = [] 
 
    #read in obs names
obs_refs,obs_locs,obs_lats,obs_lons,obs_alt,obs_number = np.loadtxt('GAW_site_indices',dtype='S20,S40,f5,f5,f5,i5',unpack=True)

#check for matching model GAW sites
def readin_gaw_sites(filename):
    reader = csv.reader(open(filename,'rb'), delimiter=',')
    for row in reader:
        new = row[:]
        try:
            locations.append(new)

        except:
            locations=[new]

    locations=np.array(locations)

    site_names = locations[:,1]

    return site_names 

model_locs = readin_gaw_sites('Simple_GAW.dat')
#strip whitespaces from model locs
model_locs = [a.replace(" ", "") for a in model_locs]

matching_sites = []

for obs in obs_locs:
	for model in model_locs:
		if obs == model:
			matching_sites = np.append(matching_sites,obs)

print matching_sites


f = 'GAW_O3_SFC_2006_2012.npy'
read = np.load(f)

times = read[0,:]
read = read[1:,:]

#iterate through sites and take LOMB
daily_mag_array = []
half_annual_mag_array = []
annual_mag_array = []

daily_phase_array = []
half_annual_phase_array = [] 
annual_phase_array = []

full_time = np.arange(0,2191,1./24)

#array containing length of months from 2006 in days

month_lengths = [31,28,31,30,31,30,31,31,30,31,30,31,  
                31,28,31,30,31,30,31,31,30,31,30,31,   
                31,29,31,30,31,30,31,31,30,31,30,31,   
                31,28,31,30,31,30,31,31,30,31,30,31,   
                31,28,31,30,31,30,31,31,30,31,30,31,    
                31,28,31,30,31,30,31,31,30,31,30,31]    

#In each iteration, limit time to the month presently in

n_months = 72

counter = 0
start_time = 0  #in mins
end_time = month_lengths[0]

for mn in range(n_months):
    #counter+=1
    #print counter
    print mn
    if mn != 0:
        start_time = end_time
        end_time = end_time + (month_lengths[mn])
    
    daily_mag_array = []
    daily_phase_array = []
    
    
    for siten in range(len(obs_refs)):
        #print siten
        vals = read[siten,:]
    
    #cut vals to monthly chunk
        chunk_test = (times >= start_time) & (times < end_time) 
        chunk_times = times[chunk_test]
        vals = vals[chunk_test]
    
        valid = vals >0
        vals = vals[valid]
        #print 'vals len =', len(vals)
        valid_times = chunk_times[valid]	

        if len(valid_times) > 0:
            window = np.hamming(len(vals))
            mean = np.mean(vals)
            vals = vals - mean
            vals = vals*window

            NOUT = 0.5*4*1*len(vals)
            NOUT = int(NOUT)

            #print NOUT
            fa, fb, mag, ph= lomb_phase.lomb(valid_times,vals,NOUT)
            periods = 1./fa
            amp_corr = 1./(sum(window)/len(window))
            mag = mag * amp_corr
	
	        #calculations for mags and phases of key periods

            closest_daily_period_index = min(range(len(periods)), key=lambda i: abs(periods[i]-daily_period))

            #print closest_daily_period_index
            daily_mag = mag[closest_daily_period_index]
            daily_phase = ph[closest_daily_period_index]
        
        else:
            daily_mag = float('NaN')
            daily_phase = float('NaN')

        daily_mag_array = np.append(daily_mag_array,daily_mag)
        daily_phase_array = np.append(daily_phase_array,daily_phase)
    
    
     
    #save fft magitude
    np.save('serial_obs_magnitudes/GAW_obs_daily_magnitudes_month%i'%(mn+1),daily_mag_array)

    #save fft phase
    np.save('serial_obs_phases/GAW_obs_daily_phases_month%i'%(mn+1),daily_phase_array)	

end = datetime.datetime.now() 
diff = end - start
print 'time = ', diff.seconds

