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

def run_LSP(valid_times,vals,start_point,end_point, x, max_size, site_ref,s_time, time_diff):
    
    vals = vals[start_point:end_point]
    valid_times = valid_times[start_point:end_point]
    
    #test if sufficient amount of data is present
    limit = (max_size/100) * 80
    if len(valid_times) >= limit:
    
        #change times to go from 0 - for accurate phase
        first_time = s_time
        valid_times = valid_times - s_time
        
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

        #print closest_daily_period_index
        daily_mag = mag[closest_daily_period_index]
        daily_phase = ph[closest_daily_period_index]
            
        #correct phase from UTC to solar time
        
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

    return (daily_mag,daily_phase, site_ref)


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

counter = 0
for obs in obs_locs:
    for model in model_locs:
        if obs == model:
            matching_sites = np.append(matching_sites,obs)
            matching_sites_n = np.append(matching_sites_n,counter)
    counter+=1
    
print matching_sites
print len(matching_sites)
print matching_sites_n

f = 'GAW_O3_SFC_2006_2012.npy'
read = np.load(f)
big_times = read[0,:]
read = read[1:,:]

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
    
    max_size = (end_time-start_time)*24
    
    big_valid_times = []
    big_vals = []
    big_start_point = []
    big_end_point = []
    big_refs = []
    time_diffs = []
    
    site_counter = 0
    for siten in matching_sites_n:
        print siten
    
        site_ref = obs_refs[siten]
        site_lon = obs_lons[siten]
        tz = timezones[siten]
        vals = read[siten,:]
        
    #cut vals to monthly chunk
        chunk_test = (big_times >= start_time) & (big_times < end_time) 
        chunk_times = big_times[chunk_test]
        
        
        vals = vals[chunk_test]
        
        valid = vals > 0
        vals = vals[valid]
        chunk_times = chunk_times[valid]
    	
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
            site_counter+=1
        else:
            start_point=end_point
            end_point+=array_len
        	 
        big_refs = np.append(big_refs,site_ref)
        big_valid_times	= np.append(big_valid_times,chunk_times)
        big_vals = np.append(big_vals,vals)
        big_start_point = np.append(big_start_point,start_point)
        big_end_point = np.append(big_end_point,end_point)
        time_diffs = np.append(time_diffs,time_diff)
        

    if __name__ == '__main__':
        pool = multiprocessing.Pool(processes=20)
        results = [pool.apply_async(run_LSP, (big_valid_times,big_vals,big_start_point[x],big_end_point[x], x, max_size,big_refs[x],start_time, time_diffs[x])) for x in range(len(matching_sites_n))]
        big_array = [r.get() for r in results]

    big_array = np.array(big_array)

    site_refs = big_array[:,2]
    
    daily_mag_array = big_array[:,0]
    daily_mag_array = np.vstack((daily_mag_array,site_refs))
    
    daily_phase_array = big_array[:,1]
    daily_phase_array = np.vstack((daily_phase_array,site_refs))
    
    pool.terminate()
     
    #save fft magitude
    np.save('obs_magnitudes/GAW_obs_daily_magnitudes_month%i'%(mn+1),daily_mag_array)

    #save fft phase
    np.save('obs_phases/GAW_obs_daily_phases_month%i'%(mn+1),daily_phase_array)	

end = datetime.datetime.now() 
diff = end - start
print 'time = ', diff.seconds

