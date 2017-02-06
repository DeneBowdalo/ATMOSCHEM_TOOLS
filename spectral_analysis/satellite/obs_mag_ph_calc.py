import numpy as np
import matplotlib.pyplot as plt
import modules
import scipy.stats
import scipy.signal
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
import numpy.fft
import glob
import sys
import multiprocessing
import operator
import collections
import re

all_data = np.load('OMI_O3_trop_ave_2x2.5.npy')
#print np.shape(all_data)

#lat_arr = np.arange(120)
#lon_arr = np.arange(288)

lat_arr = np.arange(60)
lon_arr = np.arange(144)

y = []
#reshape into stack of gridboxes
for b in lat_arr:
    print b
    for c in lon_arr:
        data = all_data[:,b,c]
	try:
            y = np.vstack((y,data))
	except:
	    y = data

#y = np.load('OMI_O3_trop_ave_gridstack.npy')
np.save('OMI_O3_trop_ave_2x2.5_gridstack',y)
#print np.shape(y)

#remove first 3 months (oct 2004 - jan 2005)
#y = y[:,3:]

#set monthly time data 
val =365.25/12
x = np.arange(len(y[0,:]))
x = x*val

def main_arg(x,y,grid_count):  
    print grid_count
    ofac = 4

    #average dt of entire time series
    diffs = [x[i+1]-x[i] for i in range(len(x)-1)]  
    avgdt = np.average(diffs)

    #make time start from 0    
    x_from0 = modules.phase_start_correct(x)

    periods,mag,ph,fr,fi,amp_corr = modules.take_lomb(x_from0,y,ofac,avgdt)

    #get mean of values
    mean_array = np.average(y)

    #correct magnitude and phase for spectral leakage 
    zoomfact = 1000
    half_annual_mag,half_annual_phase = modules.periodic_interp(fr,fi,zoomfact,periods,365.25/2.,len(y),amp_corr)
    annual_mag,annual_phase = modules.periodic_interp(fr,fi,zoomfact,periods,365.25,len(y),amp_corr)

    #correct for phase shift as data starts in Oct 2004 
    n_off = 273.25
    
    if n_off > 365.25/2:
        n_off = n_off-(365.25/2)
    offset = ((np.pi*2)/(365.25/2))*n_off
    half_annual_phase = half_annual_phase + offset
    if half_annual_phase > np.pi:
        half_annual_phase = -np.pi + (half_annual_phase - np.pi)

    n_off = 273.25
    offset = ((np.pi*2)/(365.25))*n_off
    annual_phase = annual_phase + offset
    if annual_phase > np.pi:
        annual_phase = -np.pi + (annual_phase - np.pi)

    #convert phase to time
    half_annual_phase = modules.convert_phase_units_actual_single(half_annual_phase,6)
    annual_phase = modules.convert_phase_units_actual_single(annual_phase,12)
    
    #np.save('mags_phases/mag_spectrums/%i'%(grid_count),mag)
    #np.save('mags_phases/phase_spectrums/%i'%(grid_count),ph)
    #np.save('mags_phases/periods',periods)
    return (half_annual_mag,half_annual_phase,annual_mag,annual_phase,mean_array)
    
if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=8)
    results = [pool.apply_async(main_arg, (x,y[i],i)) for i in range(len(y))]
    results = [r.get() for r in results]

results = np.array(results)
    
half_annual_mag_array = results[:,0]
half_annual_phase_array = results[:,1]
annual_mag_array = results[:,2]
annual_phase_array = results[:,3]
mean_array = results[:,4]

pool.terminate()

#save fft magitude
np.save('obs_magnitudes/obs_half_annual_magnitudes',half_annual_mag_array)
np.save('obs_magnitudes/obs_annual_magnitudes',annual_mag_array)

#save fft phase
np.save('obs_phases/obs_half_annual_phases',half_annual_phase_array)
np.save('obs_phases/obs_annual_phases',annual_phase_array)

#save average
np.save('obs_averages',mean_array)
