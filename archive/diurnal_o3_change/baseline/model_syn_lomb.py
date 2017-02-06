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
from scipy.stats import nanmean

f = '/home/db876/plotting_tools/binary_logs/GEOS_v90103_4x5_GAW_O3_logs.npy'

lon_step_time  = 24./360.

#lat lon edges for 4x5 grid 
lat_e = np.arange(-88,90,4)
lat_e = np.insert(lat_e,0,-90)
lat_e = np.append(lat_e,90)

lon_e = np.arange(-182.5,178,5)

lat_c = np.arange(-86,90,4)
lat_c = np.insert(lat_c,0,-89)
lat_c = np.append(lat_c,89)

lon_c = np.arange(-180,178,5)

#read obs data
obs_refs,obs_locs,obs_lats,obs_lons,obs_alt,timezones,obs_number = np.loadtxt('GAW_site_indices',dtype='S20,S20,f5,f5,f5,i5,i5',unpack=True)
#get site name from ref
o_test = obs_refs == 'smo'
site_name = obs_locs[o_test]
site_name = site_name[0]
tz = timezones[o_test]

#get obs_lon of site
site_lon = obs_lons[o_test]
site_num = obs_number[o_test]
site_num = site_num[0]

#obs_lon = obs_lon[0]
read = np.load(f)
read = read[1:,:]

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


for obs in obs_locs:
    counter = 0
    for model in model_locs:
        if obs == model:
            matching_sites = np.append(matching_sites,obs)
            matching_sites_n = np.append(matching_sites_n,counter)
        counter+=1

matching_sites_n = np.insert(matching_sites_n,8,99999)

print matching_sites
print len(matching_sites_n)

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


daily_period = 1

#hourly_step = 1./24

#a = np.arange(0,25,hourly_step)
#b = np.arange(9,78,hourly_step)


#pi2 = np.pi*2

#daily_amplitude = 2


#cos_waveform1 = daily_amplitude*(np.cos((pi2*b/5)-(0)))
#vals = 50 + (cos_waveform1)

month_lengths = [31,28,31,30,31,30,31,31,30,31,30,31,  
                31,28,31,30,31,30,31,31,30,31,30,31,   
                31,29,31,30,31,30,31,31,30,31,30,31,   
                31,28,31,30,31,30,31,31,30,31,30,31,   
                31,28,31,30,31,30,31,31,30,31,30,31,    
                31,28,31,30,31,30,31,31,30,31,30,31]  

start = 0
end = 31*24

counter = 1

for n in month_lengths:

    
    ave_day = []
    print site_num
    print (matching_sites_n[site_num-1]) +1
    valid = labels == (matching_sites_n[site_num-1]) +1
    vals =  all_vals[valid]
    vals = vals*1e9
    
    vals = vals[start:end]

    other_vals = vals[:]
    
    #print np.min(vals)
    
    for i in range(len(vals)):
        if vals[i] <= 0:
            other_vals[i] = np.nan
    
    b = np.arange(0,n,1./24)
    
    #convert site_lon to 0 to 360 degs
    if site_lon < 0:
    	site_lon = 360-np.abs(site_lon)
    	
    #transform from UTC time to sun time 
    sun_time = lon_step_time*site_lon
    time_diff = sun_time - 0
    if time_diff > 12:
    	time_diff = time_diff-24
    
    max_size = n * 24
    
    limit = (max_size/100) * 80
    print len(vals)
    if len(b) >= limit:
    
        for i in range(24):
            a_day = other_vals[i::24]
            ave_day = np.append(ave_day,nanmean(a_day))
            
        if nanmean(ave_day) > 0:
            fig=plt.figure(figsize=(20,12))
            fig.patch.set_facecolor('white')
            xi = np.arange(0,24,1)
            #print 'ave_day= ', ave_day
            #print len(ave_day)
            #plt.plot(range(len(other_vals)),other_vals)	
            plt.plot(xi,ave_day)
            print len(xi)
        #plt.loglog(periods, mag, color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='Data')	
           
        
        
        
        
        window = np.hamming(len(vals))
        mean = np.mean(vals)
        vals = vals - mean
        vals = vals*window


        NOUT = 0.5*4*1*len(vals)
        NOUT = int(NOUT)

        fa, fb, mag, ph= lomb_phase.lomb(b,vals,NOUT)
        periods = 1./fa
        amp_corr = 1./(sum(window)/len(window))
        mag = mag * amp_corr
	

        closest_daily_period_index = min(range(len(periods)), key=lambda i: abs(periods[i]-daily_period))
      
        daily_mag = mag[closest_daily_period_index]
        daily_phase = ph[closest_daily_period_index]
        
        print np.min(ph)
        print np.max(ph)
        #print vals
        print 'daily phase = ', daily_phase
        
        print 'daily phase = ', modules.convert_phase_units_actual_single(daily_phase,24)
        
        #daily_phase = modules.correct_select_daily_phase_actual_single(daily_phase,lon_c,lat_c,lon_e,lat_e, site_lon)
        
        #daily_phase = modules.convert_phase_units_actual_single(daily_phase,24)
        
    
        print 'sun time =', sun_time
        print 'time diff = ', time_diff
        #adjust daily_phase to solar time
        correction_f = ((2*np.pi)/24) * time_diff
        daily_phase = daily_phase + (correction_f)
        if daily_phase < -np.pi:
            diff = np.abs(daily_phase) - np.pi
            daily_phase = np.pi - diff 
        if daily_phase > np.pi:
            diff = np.abs(daily_phase) - np.pi
            daily_phase = -np.pi + diff
        
        print 'daily phase = ', modules.convert_phase_units_actual_single(daily_phase,24)
        
        #daily_phase = -np.pi
        #daily_phase = modules.correct_select_daily_phase_actual_single(daily_phase,lon_c,lat_c,lon_e,lat_e, site_lon)
        daily_phase = modules.convert_phase_units_actual_single(daily_phase,24)
        
        print 'daily phase =', daily_phase
        
        
        print 'daily mag = ', daily_mag
        plt.show()
    else:
        daily_mag = np.nan
        daily_phase = np.nan
        print 'skip'
        
    start=end
    end+=((month_lengths[counter])*24)
    counter+=1
    
    

    #print 'daily magnitude', daily_mag
    