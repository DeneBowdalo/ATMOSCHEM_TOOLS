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
import redfit
import numpy.fft


#set parameters
mctest = False
nsim= 100
n50 = 1
ofac = 1

start = datetime.datetime.now()
daily_period = 1
half_annual_period = 365.25/2
annual_period = 365.25

lon_step_time  = 24./360.

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-3]
years = paths[-2]
start_year = years[:4]
end_year = years[5:9]
model_path = paths[-1] 

data_split = model_path.split('_SFC') 
print data_split
model = data_split[0]
#model = model[4:]

only_model = False
model_version = False
model_version_grid = False
model_version_grid_met = False

other = data_split[1]
other_split = other.split('_')
other_split.remove('')
if other_split == []:
    print 'only 1 version of model'
    only_model = True
elif len(other_split) == 1:
    version = other_split[0]
    model_version = True
elif len(other_split) == 2:
    version = other_split[0]
    grid = other_split[1]
    model_version_grid = True
elif len(other_split) == 3:
    version = other_split[0]
    grid = other_split[1]    
    met = other_split[2]
    model_version_grid_met = True

print start_year, end_year
print model 


print '\nSpecies is %s\n'%(species)


def run_LSP(mod_data,x):

    lat_i = lat_indices[x]
    lon_i = lon_indices[x]

    print lat_i,lon_i

    current_lat = lat_c[lat_i]
    current_lon = lon_c[lon_i]

    waveform = mod_data
    
    waveform_ave = np.average(waveform)
    
    model_date_val = np.copy(model_date)
    model_time_val = np.copy(model_time)
    
    time = modules.date_process(model_date_val,model_time_val,start_year)
    
    if (species.lower() != 'gmao_temp') and (species.lower() != 'gmao_psfc') and (species.lower() != 'wind_speed') and (species.lower() != 'wind_direction'):
        waveform = waveform*1e9	
 
    #check model vals are valid
    #valid = vals >= 0
    #vals = vals[valid]
    #model_time_val = model_time[valid]
    #model_date_val = model_date[valid]

    #take 8 hour average
    divisor = 8

    total_len = len(waveform)/divisor
    start = 0
    end = divisor
    ave_waveform = []
    ave_time = []
    for i in range(total_len):
        ave = np.ma.average(waveform[start:end])
        ave_time=np.append(ave_time,time[start])
        ave_waveform=np.append(ave_waveform,ave)
        start+=divisor
        end+=divisor
 
    time=np.copy(ave_time)
    waveform=np.copy(ave_waveform)

    #take lsp unwindowed of waveform

    ua_periods,ua_mag,ua_ph,ua_fr,ua_fi = modules.take_lomb_unwindowed(time,waveform,ofac,1./24)

    #take out known periodic components 1,182.625, and 365.25 a priori for more accurate red noise fit.
    closest_daily_index = min(range(len(ua_periods)), key=lambda i: abs(ua_periods[i]-1.))
    closest_ha_index = min(range(len(ua_periods)), key=lambda i: abs(ua_periods[i]-182.625))
    closest_annual_index = min(range(len(ua_periods)), key=lambda i: abs(ua_periods[i]-365.25))

    rm_indices = [closest_daily_index,closest_ha_index,closest_annual_index]

    ua_mag_c,ua_fr,ua_fi = redfit.sidelobe_percent_remove(np.copy(ua_mag),ua_fr,ua_fi,rm_indices,5.,ua_periods)
    
    #-------------------------------------------------------------------------------
    #Do IFFT of altered spectra - with significant periods removed and gaps left in real and imag components linearly interpolated.
    #altered spectra provides red noise estimation baseline

    ##use ifft to get time series back from adjusted spectra
    #complex Fourier spectrum which corresponds to the Lomb-Scargle periodogram:
    F = [0]*((len(ua_fr)*2)+1)

    #set first real value to average 
    F[0] = complex(waveform_ave*len(waveform),0)

    #Get reverse real and imaginary values
    rev_ua_fr=np.copy(ua_fr[::-1])
    rev_ua_fi=np.copy(ua_fi[::-1])

    rev_ua_fr[0] = 0
    rev_ua_fi[0] = 0

    f_index = 1

    #Fill Fourier Spectrum real and imaginary values
    for i in range(len(ua_fr)):
        F[f_index] = complex(ua_fr[i],ua_fi[i])
        f_index+=1

    for i in range(len(ua_fr)):
        F[f_index] = complex(rev_ua_fr[i],-rev_ua_fi[i])
        f_index+=1

    F = np.array(F)    

    #Take ifft and just take real values
    ifft_ua_ts = numpy.fft.ifft(F)
    ifft_ua_ts = ifft_ua_ts.astype('float64')

    ifft_ua_ts_len = (len(ifft_ua_ts)/ofac) + np.mod(len(ifft_ua_ts),ofac)

    ifft_time = time[-ifft_ua_ts_len:]
    ifft_ua_ts = ifft_ua_ts[-len(waveform):]

    ifft_time = ifft_time-ifft_time[0]

    a_periods,a_mag,corr_a_mag,a_fr,a_fi,a_red_periods,a_red_mag,a_gredth,a_fac95,a_fac99,a_fac99_9,a_faccrit,a_fac_grid,a_sig_levels,a_tau,a_corr = redfit.red_background(nsim,mctest,ifft_time,ifft_ua_ts,ofac)

    #apply lsp correction from altered spectrum to unaltered spectrum
    corr_ua_mag = ua_mag/a_corr

    #check confidence of each point on spectrum

    sigs = np.zeros(len(corr_ua_mag))

    last_ind = len(a_sig_levels)-1

    for i in range(len(a_sig_levels)-1):
        conf_low = a_gredth*a_fac_grid[i]
        conf_up = a_gredth*a_fac_grid[i+1]
    
        current_last_ind = i+1
    
        for j in range(len(corr_ua_mag)):
            if sigs[j] == 0:
                if (corr_ua_mag[j] >= conf_low[j]) and (corr_ua_mag[j] < conf_up[j]):
                    sigs[j] = a_sig_levels[i]
                elif current_last_ind == last_ind:
                    if corr_ua_mag[j] > conf_up[j]:
                       sigs[j] = a_sig_levels[i+1]
    
    #get critical significance for all points on spectrum
    crit_sig = a_gredth*a_faccrit
    
    #get 95,99 and 99.9 % chi squared significance bands for all points on spectrum
    sig_95 = a_gredth*a_fac95
    sig_99 = a_gredth*a_fac99
    sig_99_9 = a_gredth*a_fac99_9
    
    return (x,sigs,sig_95,sig_99,sig_99_9,crit_sig,a_gredth,corr_ua_mag,ua_periods,a_tau)

#read model netcdf file


if only_model == True:
    model_root_grp = Dataset('/work/home/db876/plotting_tools/model_files/%s_SURFACE_%s_%s_%s.nc'%(model,species,start_year,end_year))
elif model_version == True:
    model_root_grp = Dataset('/work/home/db876/plotting_tools/model_files/%s_SURFACE_%s_%s_%s_%s.nc'%(model,species,start_year,end_year,version))
elif model_version_grid == True:
    model_root_grp = Dataset('/work/home/db876/plotting_tools/model_files/%s_SURFACE_%s_%s_%s_%s_%s.nc'%(model,species,start_year,end_year,version,grid))
elif model_version_grid_met == True:
    model_root_grp = Dataset('/work/home/db876/plotting_tools/model_files/%s_SURFACE_%s_%s_%s_%s_%s_%s.nc'%(model,species,start_year,end_year,version,grid,met))

model_data = model_root_grp.variables[species.lower()][:]
model_date = model_root_grp.variables['date'][:]
model_time = model_root_grp.variables['time'][:]
lat_c = model_root_grp.variables['lat_centre'][:]
lon_c = model_root_grp.variables['lon_centre'][:]
lat_e = model_root_grp.variables['lat_edges'][:]
lon_e = model_root_grp.variables['lon_edges'][:]

n_boxes = len(lat_c)*len(lon_c)
#return gridbox number required from lat lon of obs

lat_indices = []
lon_indices = []
linear_data = []    

lat_i = 0
lon_i = 0

#n_boxes = 10
#lat_c = lat_c[:2]
#lon_c = lon_c[:5]
#lat_e = lat_c[:2]
#lon_e = lon_c[:5]

for siten in range(n_boxes):
    linear_data.append(model_data[:,lat_i,lon_i])
   
    lat_indices.append(lat_i)
    lon_indices.append(lon_i)
   
    if lon_i == (len(lon_c)-1):
        lat_i+=1
        lon_i=0
    else:
        lon_i+=1

if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=16)
    results = [pool.apply_async(run_LSP, (linear_data[x],x)) for x in range(n_boxes)]
    big_array = [r.get() for r in results]
    pool.terminate()

indices_array = []
chi_sig_array = []
chi95_sig_array = []
chi99_sig_array = []
chi99_9_sig_array = []
critical_sig_array = []
redfit_array = []
mag_array = []
period_array = []
tau_array = []

for i in range(len(big_array)):
    cut = big_array[i]
    
    indices_array.append(cut[0])
    chi_sig_array.append(cut[1])
    chi95_sig_array.append(cut[2])
    chi99_sig_array.append(cut[3])
    chi99_9_sig_array.append(cut[4])
    critical_sig_array.append(cut[5])
    redfit_array.append(cut[6])
    mag_array.append(cut[7])
    period_array.append(cut[8])
    tau_array.append(cut[9])

chi_sig_array = np.array(chi_sig_array)
chi95_sig_array = np.array(chi95_sig_array)
chi99_sig_array = np.array(chi99_sig_array)
chi99_9_sig_array = np.array(chi99_9_sig_array)
critical_sig_array = np.array(critical_sig_array)
redfit_array = np.array(redfit_array)
mag_array = np.array(mag_array)
period_array = np.array(period_array)
tau_array = np.array(tau_array)

#sort arrays by indices array for sanity
chi_sig_array = chi_sig_array[indices_array]
chi95_sig_array = chi95_sig_array[indices_array]
chi99_sig_array = chi99_sig_array[indices_array]
chi99_9_sig_array = chi99_9_sig_array[indices_array]
critical_sig_array = critical_sig_array[indices_array]
redfit_array = redfit_array[indices_array]
mag_array = mag_array[indices_array]
period_array = period_array[indices_array]
tau_array = tau_array[indices_array]

#save out spectra to netcdf
root_grp_spec = Dataset('model_spectra_significance.nc', 'w')
root_grp_spec.description = 'Significance for spectra from Lomb-Scargle Periodogram analysis of Model - Program written by Dene Bowdalo'
root_grp_spec.createDimension('site_n', None)
root_grp_spec.createDimension('lat_centre', len(lat_c))
root_grp_spec.createDimension('lon_centre', len(lon_c))
root_grp_spec.createDimension('lat_edges', len(lat_e))
root_grp_spec.createDimension('lon_edges', len(lon_e))
root_grp_spec.createDimension('grid_s', 1)
root_grp_spec.createDimension('n_boxes', n_boxes)

chi_sig_array = np.reshape(chi_sig_array,(len(lat_c),len(lon_c),-1))
chi95_sig_array = np.reshape(chi95_sig_array,(len(lat_c),len(lon_c),-1))
chi99_sig_array = np.reshape(chi99_sig_array,(len(lat_c),len(lon_c),-1))
chi99_9_sig_array = np.reshape(chi99_9_sig_array,(len(lat_c),len(lon_c),-1))
critical_sig_array = np.reshape(critical_sig_array,(len(lat_c),len(lon_c),-1))
redfit_array = np.reshape(redfit_array,(len(lat_c),len(lon_c),-1))
mag_array = np.reshape(mag_array,(len(lat_c),len(lon_c),-1))
tau_array = np.reshape(tau_array,(len(lat_c),len(lon_c)))

chi_significance = root_grp_spec.createVariable('chi_significance_points', 'f8', ('lat_centre','lon_centre','site_n'))
chi_significance_95 = root_grp_spec.createVariable('chi_significance_95', 'f8', ('lat_centre','lon_centre','site_n'))
chi_significance_99 = root_grp_spec.createVariable('chi_significance_99', 'f8', ('lat_centre','lon_centre','site_n'))
chi_significance_99_9 = root_grp_spec.createVariable('chi_significance_99_9', 'f8', ('lat_centre','lon_centre','site_n'))
critical_significance = root_grp_spec.createVariable('critical_significance', 'f8', ('lat_centre','lon_centre','site_n'))
red_fit = root_grp_spec.createVariable('red_noise_fit', 'f8', ('lat_centre','lon_centre','site_n'))
amplitude = root_grp_spec.createVariable('amplitude', 'f8', ('lat_centre','lon_centre','site_n'))
period = root_grp_spec.createVariable('period', 'f8', ('site_n'))
tau = root_grp_spec.createVariable('tau', 'f8', ('lat_centre','lon_centre'))
lat_centre = root_grp_spec.createVariable('lat_centre', 'f8', ('lat_centre',))
lon_centre = root_grp_spec.createVariable('lon_centre', 'f8', ('lon_centre',))
lat_edge = root_grp_spec.createVariable('lat_edges', 'f8', ('lat_edges',))
lon_edge = root_grp_spec.createVariable('lon_edges', 'f8', ('lon_edges',))   
    
chi_significance[:] = chi_sig_array
chi_significance_95[:] = chi95_sig_array
chi_significance_99[:] = chi99_sig_array
chi_significance_99_9[:] = chi99_9_sig_array
critical_significance[:] = critical_sig_array
red_fit[:] = redfit_array
amplitude[:] = mag_array
period[:] = period_array[0]
tau[:] = tau_array
lat_centre[:] = lat_c
lon_centre[:] = lon_c
lat_edge[:] = lat_e     
lon_edge[:] = lon_e       

#plt.loglog(period_array[0],mag_array[0,0,:])
#plt.loglog(period_array[0],red_fit[0,0,:])
#plt.loglog(period_array[0],chi_significance_99_9[0,0,:])
#plt.show()

root_grp_spec.close()

end = datetime.datetime.now() 
diff = end - start
print 'time = ', diff.seconds
