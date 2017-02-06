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
from netCDF4 import Dataset
from scipy import stats

present_dir = os.getcwd()

species = 'O3'
start_year = int(present_dir[-2][:4])
end_year = int(present_dir[-2][-4:])
start_year_file = int(present_dir[-3][:4])
end_year_file = int(present_dir-3][-4:])

n_years = end_year-start_year

start = datetime.datetime.now()
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.

def run_LSP(obs_time,site_lon,x):

    data_valid = True

    #print 'lev %s'%(x)
    
    full_times_year = np.arange(0,365,1.)
    
    #check obs vals are valid
    valid = vals >= 0
    vals = vals[valid]
    valid_times = obs_time[valid]

    #if length of vals is zero then class as invalid immediately
    if len(vals) == 0:
        data_valid = False
    else:
        #test if there if data is valid to process at each height for each site
        #data should not have gaps > 1 year or 
        time_gaps = np.diff(valid_times)                
        inv_count = 0
        max_count = round(n_years/2.)
        for i in time_gaps:
            if i > 90:
                inv_count+=1
                if inv_count >= max_count:
                    data_valid = False
                    print 'Persisent Data gap > 3 months'
                    break
            if i > 365:
                    data_valid = False
                    print 'Data gap > 1 Year'
                    break

    if data_valid == True:
    
        #convert site_lon to 0 to 360 degs
        if site_lon < 0:
            site_lon = 360-np.abs(site_lon)

        #make time start from 0    
        valid_times_from0 = modules.phase_start_correct(valid_times)

        periodic_periods = [365.25/4.,365.25/3.,365.25/2.,365.25]
        periods,mag,ph,fr,fi = modules.take_lomb_spec(valid_times_from0,vals,w=True,key_periods=periodic_periods)

        #get mean of values
        mean_array = np.average(vals)
    
        #correct all phases for start point (not actually being from 0 - just corrected to be)
        ph = modules.phase_start_point_correct_all(periodic_periods,ph,valid_times)

        key_seasonal_periods = [365.25/4.,365.25/3.,365.25/2.,365.25]

        seasonal_mags = mag[:]
        seasonal_phs = ph[:]

        seasonal_h3_mag = mag[0]
        seasonal_h2_mag = mag[1]
        seasonal_h1_mag = mag[2]
        annual_mag = mag[3]
        seasonal_h3_ph = ph[0]
        seasonal_h2_ph = ph[1]
        seasonal_h1_ph = ph[2]
        annual_ph = ph[3]

        #convolve annual cycle and harmonics to seasonal waveform for 1 year
        seasonal_mag,seasonal_min_ph,seasonal_max_ph,seasonal_waveform,seasonal_ff = modules.period_convolution(key_seasonal_periods,full_times_year,seasonal_mags,seasonal_phs,mean_array)
        
        #convert phase to time
        seasonal_h3_ph = modules.convert_phase_units_actual_single(seasonal_h3_ph,3.)
        seasonal_h2_ph = modules.convert_phase_units_actual_single(seasonal_h2_ph,4.)
        seasonal_h1_ph = modules.convert_phase_units_actual_single(seasonal_h1_ph,6.)
        annual_ph = modules.convert_phase_units_actual_single(annual_ph,12.)
        seasonal_min_ph = modules.convert_phase_units_actual_single(seasonal_min_ph,12.)
        seasonal_max_ph = modules.convert_phase_units_actual_single(seasonal_max_ph,12.)

    else:
        seasonal_h3_mag = -99999
        seasonal_h2_mag = -99999
        seasonal_h1_mag = -99999
        annual_mag = -99999
        seasonal_mag = -99999
        
        seasonal_h3_ph = -99999
        seasonal_h2_ph = -99999
        seasonal_h1_ph = -99999
        annual_ph = -99999
        seasonal_max_ph = -99999
        seasonal_min_ph = -99999

        seasonal_waveform = np.array([-99999]*len(full_times_year))

        mean_array = -99999
    
    return x,seasonal_h3_mag,seasonal_h3_ph,seasonal_h2_mag,seasonal_h2_ph,seasonal_h1_mag,seasonal_h1_ph,annual_mag,annual_ph,seasonal_mag,seasonal_max_ph,seasonal_min_ph,seasonal_waveform,mean_array


obs_file = '/work/home/db876/observations/ozonesonde/process/%s_RADIOSONDES_%s_%s.nc'%(species,start_year_file,end_year_file)
obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,gap_inds = modules.read_obs_all(obs_file,species,start_year,end_year)

#read in obs names
root_grp = Dataset(obs_file)
    
pressure_centres = np.empty((len(obs_refs),47))

for i in range(len(obs_refs)):
    ref_grp = root_grp.groups[obs_refs[i]]
    pressure_centres[i,:] = ref_grp.variables['pressure_centres'][:]

n_levs = 47

#setup netcdf
root_grp_out = Dataset('obs_sig_periods.nc', 'w')
root_grp_out.description = 'Amplitudes and Phases for key periods from Lomb-Scargle Periodogram analysis of Ozonesonde Observations - Program written by Dene Bowdalo'
root_grp_out.createDimension('pressures', n_levs)
root_grp_out.createDimension('refs', len(obs_refs))
root_grp_out.createDimension('sw', 365)

sha3 = root_grp_out.createVariable('seasonal_harmonic3_amplitude','f8',('refs','pressures'))    
sha2 = root_grp_out.createVariable('seasonal_harmonic2_amplitude','f8',('refs','pressures')) 
sha1 = root_grp_out.createVariable('seasonal_harmonic1_amplitude','f8',('refs','pressures'))
aa = root_grp_out.createVariable('annual_amplitude','f8',('refs','pressures'))
sa = root_grp_out.createVariable('seasonal_amplitude','f8',('refs','pressures'))
shp3 = root_grp_out.createVariable('seasonal_harmonic3_phase','f8',('refs','pressures'))    
shp2 = root_grp_out.createVariable('seasonal_harmonic2_phase','f8',('refs','pressures')) 
shp1 = root_grp_out.createVariable('seasonal_harmonic1_phase','f8',('refs','pressures'))
ap = root_grp_out.createVariable('annual_phase','f8',('refs','pressures'))
smaxp = root_grp_out.createVariable('seasonal_max_phase','f8',('refs','pressures'))
sminp = root_grp_out.createVariable('seasonal_min_phase','f8',('refs','pressures'))
sw = root_grp_out.createVariable('seasonal_waveform','f8',('refs','sw','pressures'))
ave = root_grp_out.createVariable('average','f8',('refs','pressures'))  
refs = root_grp_out.createVariable('refs',str,('refs'))
p_c = root_grp_out.createVariable('pressure_centres','f8',('refs','pressures'))
lat = root_grp_out.createVariable('latitude','f8',('refs',))
lon = root_grp_out.createVariable('longitude','f8',('refs',))
cnt = root_grp_out.createVariable('countries',str,('refs',))

refs[:] = np.array(obs_refs,dtype=object) 
p_c[:] = pressure_centres
lat[:] = lats
lon[:] = lons
cnt[:] = np.array(countries,dtype=object)

for j in range(len(obs_refs)):
    ref = obs_refs[j]
    print ref

    linear_data = []    

    site_group = root_grp.groups[ref] 
    data = site_group.variables['o3'][:]
    time = site_group.variables['time'][:] 
    time = time - time[0]
    
    for i in range(n_levs):
        linear_data.append(data[:,i])  

    obs_lon = site_group.longitude
    
    if __name__ == '__main__':
        pool = multiprocessing.Pool(processes=16)
        results = [pool.apply_async(run_LSP, (time,linear_data[x],obs_lon,x)) for x in range(n_levs)]
        big_array = [r.get() for r in results]

        pool.terminate()

    indices_array = []
    seasonal_h3_mag_array = []
    seasonal_h2_mag_array = []
    seasonal_h1_mag_array = []
    annual_mag_array = []
    seasonal_mag_array = []
    seasonal_h3_ph_array = []
    seasonal_h2_ph_array = []
    seasonal_h1_ph_array = []
    annual_ph_array = []
    seasonal_ph_max_array = []
    seasonal_ph_min_array = []
    seasonal_waveform = []
    mean_array = []

    for i in range(len(big_array)):
        cut = big_array[i]
        
        indices_array.append(cut[0])
        seasonal_h3_mag_array.append(cut[1])
        seasonal_h2_mag_array.append(cut[2])
        seasonal_h1_mag_array.append(cut[3])
        annual_mag_array.append(cut[4])
        seasonal_mag_array.append(cut[5])
        seasonal_h3_ph_array.append(cut[6])
        seasonal_h2_ph_array.append(cut[7])
        seasonal_h1_ph_array.append(cut[8])
        annual_ph_array.append(cut[9])
        seasonal_ph_max_array.append(cut[10])
        seasonal_ph_min_array.append(cut[11])
        seasonal_waveform.append(cut[12])
        mean_array.append(cut[13])
    
    seasonal_h3_mag_array = np.array(seasonal_h3_mag_array)
    seasonal_h2_mag_array = np.array(seasonal_h2_mag_array)
    seasonal_h1_mag_array = np.array(seasonal_h1_mag_array)
    annual_mag_array = np.array(annual_mag_array)
    seasonal_h3_ph_array = np.array(seasonal_h3_ph_array)
    seasonal_h2_ph_array = np.array(seasonal_h2_ph_array)
    seasonal_h1_ph_array = np.array(seasonal_h1_ph_array)
    annual_ph_array = np.array(annual_ph_array)
    seasonal_ph_max_array = np.array(seasonal_ph_max_array)
    seasonal_ph_min_array = np.array(seasonal_ph_min_array)
    seasonal_waveform  = np.array(seasonal_waveform)
    mean_array = np.array(mean_array)

    #sort arrays by indices array for sanity
    seasonal_h3_mag_array = seasonal_h3_mag_array[indices_array]
    seasonal_h2_mag_array = seasonal_h2_mag_array[indices_array]
    seasonal_h1_mag_array = seasonal_h1_mag_array[indices_array]
    annual_mag_array = annual_mag_array[indices_array]
    seasonal_h3_ph_array = seasonal_h3_ph_array[indices_array]
    seasonal_h2_ph_array = seasonal_h2_ph_array[indices_array]
    seasonal_h1_ph_array = seasonal_h1_ph_array[indices_array]
    annual_ph_array = annual_ph_array[indices_array]
    seasonal_ph_max_array = seasonal_ph_max_array[indices_array]
    seasonal_ph_min_array = seasonal_ph_min_array[indices_array]
    #seasonal_waveform = seasonal_waveform[indices_array]
    mean_array = mean_array[indices_array]
 
    sha3[j,:] = seasonal_h3_mag_array
    sha2[j,:] = seasonal_h2_mag_array
    sha1[j,:] = seasonal_h1_mag_array
    aa[j,:] = annual_mag_array
    sa[j,:] = seasonal_mag_array
    shp3[j,:] = seasonal_h3_ph_array
    shp2[j,:] = seasonal_h2_ph_array
    shp1[j,:] = seasonal_h1_ph_array
    ap[j,:] = annual_ph_array
    smaxp[j,:] = seasonal_ph_max_array
    sminp[j,:] = seasonal_ph_min_array
    
    for i in range(len(seasonal_waveform)):
        sw[j,:,i] = seasonal_waveform[i]
    
    ave[j,:] = mean_array

root_grp_out.close()

end = datetime.datetime.now() 
diff = end - start
print 'time = ', diff.seconds

