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
import modules

start = datetime.datetime.now()

daily_period = 1

lon_step_time  = 24./360.

def run_LSP(valid_times,vals,start_point,end_point, x, time_diff):
    
    vals = vals[start_point:end_point]
    valid_times = valid_times[start_point:end_point]
    
    if len(valid_times) > 0:
    
        #change times to go from 0 - for accurate phase
        first_time = valid_times[0]
        valid_times = valid_times - first_time
    
        #window
        window = np.hamming(len(vals))
        mean = np.mean(vals)
        vals = vals - mean
        vals = vals*window

        NOUT = 0.5*4*1*len(vals)
        NOUT = int(NOUT)

        #take lomb
        fa, fb, mag, ph= lomb_phase.lomb(valid_times,vals,NOUT)
        periods = 1./fa
        amp_corr = 1./(sum(window)/len(window))
        mag = mag * amp_corr
	
	    #calculations for mags and phases of key periods

        closest_daily_period_index = min(range(len(periods)), key=lambda i: abs(periods[i]-daily_period))

        daily_mag = mag[closest_daily_period_index]
        daily_phase = ph[closest_daily_period_index]
        
       #adjust daily_phase to solar time
        correction_f = ((2*np.pi)/24) * time_diff
        daily_phase = daily_phase + (correction_f)
        if daily_phase < -np.pi:
            diff = np.abs(daily_phase) - np.pi
            daily_phase = np.pi - diff 
        if daily_phase > np.pi:
            diff = np.abs(daily_phase) - np.pi
            daily_phase = -np.pi + diff
            
        #convert phase to hours    
        daily_phase = modules.convert_phase_units_actual_single(daily_phase,24)

    else:
        print 'nan'
        daily_mag = float('NaN')
        daily_phase = float('NaN')

    return (daily_mag,daily_phase, x)


 #read in obs names
obs_refs,obs_locs,obs_lats,obs_lons,obs_alt,timezones,obs_number = np.loadtxt('GAW_site_indices',dtype='S20,S40,f5,f5,f5,i5,i5',unpack=True)

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
matching_sites_n = []
site_refs = []
site_lons = []

counter2 = 0
for obs in obs_locs:
    counter = 0
    for model in model_locs:
        if obs == model:
            matching_sites = np.append(matching_sites,obs)
            matching_sites_n = np.append(matching_sites_n,counter)
            site_refs = np.append(site_refs,obs_refs[counter2])
            site_lons = np.append(site_lons,obs_lons[counter2])
        counter+=1
    counter2+=1
        
f = '/home/db876/plotting_tools/binary_logs/drydepoff_gaw_logs_O3.npy'

read = np.load(f)
read = read[1:,:]

labels = read[:,0]
labels = labels.astype(int)

valid = labels == 1
dates = read[valid,1]
dates = dates.astype(int)
hours = read[valid,2]
hours = hours.astype(int)

all_vals =  read[:,3]
all_vals = all_vals.astype(float)

big_times = modules.date_process(dates,hours)
big_times = np.array(big_times)

#iterate through sites and take LOMB
daily_mag_array = []
daily_phase_array = []

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
start_time = 0  #in days
end_time = month_lengths[0]

for mn in range(n_months):
    print mn
    if mn != 0:
        start_time = end_time
        end_time = end_time + (month_lengths[mn])
    
    big_valid_times = []
    big_vals = []
    big_start_point = []
    big_end_point = []
    time_diffs = []
    
    site_counter = 0
    for siten in matching_sites_n:
        print siten
        
        valid = labels == siten+1
        vals =  all_vals[valid]
        vals = vals*1e9
        
        site_ref = site_refs[site_counter]
        site_lon = site_lons[site_counter]
        
    #cut vals to monthly chunk
        chunk_test = (big_times >= start_time) & (big_times < end_time) 
        chunk_times = big_times[chunk_test]
        vals = vals[chunk_test]
        
        
        #convert site_lon to 0 to 360 degs
    	if site_lon < 0:
    	    site_lon = 360-np.abs(site_lon)
    	
    	#transform from UTC time to solar time 
    	sun_time = lon_step_time*site_lon
    	time_diff = sun_time - 0
    	if time_diff > 12:
    	    time_diff = time_diff-24
    	         
    	#set start and point of multiple site arrays     
        array_len = len(vals)
        if site_counter == 0:
            start_point = 0
            end_point = array_len
        else:
            start_point=end_point
            end_point+=array_len
        	
        big_valid_times	= np.append(big_valid_times,chunk_times)
        big_vals = np.append(big_vals,vals)
        big_start_point = np.append(big_start_point,start_point)
        big_end_point = np.append(big_end_point,end_point)
        time_diffs = np.append(time_diffs,time_diff)
        

        site_counter+=1

    if __name__ == '__main__':
        pool = multiprocessing.Pool(processes=20)
        results = [pool.apply_async(run_LSP, (big_valid_times,big_vals,big_start_point[x],big_end_point[x], x, time_diffs[x])) for x in range(len(matching_sites_n))]
        big_array = [r.get() for r in results]
        #pool.close()
        #pool.join()

    big_array = np.array(big_array)

    daily_mag_array = big_array[:,0]
    daily_phase_array = big_array[:,1]
       
    #print 'daily phase array = ',daily_phase_array
    
    pool.terminate()
     
    #save fft magitude
    np.save('model_magnitudes/GAW_model_daily_magnitudes_month%i'%(mn+1),daily_mag_array)

    #save fft phase
    np.save('model_phases/GAW_model_daily_phases_month%i'%(mn+1),daily_phase_array)	

end = datetime.datetime.now() 
diff = end - start
print 'time = ', diff.seconds

