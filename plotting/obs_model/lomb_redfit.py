import numpy as np
import matplotlib.pyplot as plt
import glob
import lomb_phase
import modules
import datetime
from netCDF4 import Dataset
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
from scipy.stats import norm
import scipy.signal
import redfit
import os
import numpy.fft

#get species from current directory
present_dir = os.getcwd()
paths = present_dir.split("/")
species = paths[-2]
print '\nSpecies is %s\n'%(species)

#process model
#----------------------------------------
#read in model data

model_dict = {}
version_dict = {}
grid_dict = {}
met_dict = {}

#find valid_models
model_files = glob.glob('/work/home/db876/plotting_tools/model_files/*%s*'%(species))
all_models = []
for i in model_files:
    i = i.replace("/work/home/db876/plotting_tools/model_files/", "")
    model_split = i.split('_%s_'%(species))
    part1 = model_split[0]
    part2 = model_split[1]
    start_year = part2[:4]
    end_year = part2[5:9]
    year_range = start_year+'_'+end_year
    model = part1[:-8]
    #check for version, grid_size and metorology
    #if only 1 there is version and grid_size, if 2 there is version and grid_size, if 3 there it is version, grid_size and meterology.
    extra_data = part2[10:-3]
    extra_data_split = extra_data.split('_')
    
    if (len(extra_data_split) == 1) & (extra_data_split != ['']):
        version = extra_data_split[0]
        try:
            key = version_dict[model]
            if version not in key:
                key.append(version) 
                version_dict[model] = key
        except:
            version_dict[model] = [version]
    elif len(extra_data_split) == 2:
        version = extra_data_split[0]
        grid = extra_data_split[1]
        try:
            key = version_dict[model]
            if version not in key:
                key.append(version) 
                version_dict[model] = key
        except:
            version_dict[model] = [version]
        try:
            key = grid_dict[model]
            if grid not in key:
                key.append(grid) 
                grid_dict[model] = key
        except:
            grid_dict[model] = [grid]

    elif len(extra_data_split) == 3:
        version = extra_data_split[0]
        grid = extra_data_split[1]
        met = extra_data_split[2]
        try:
            key = version_dict[model]
            if version not in key:
                key.append(version) 
                version_dict[model] = key
        except:
            version_dict[model] = [version]
        try:
            key = grid_dict[model]
            if grid not in key:
                key.append(grid) 
                grid_dict[model] = key
        except:
            grid_dict[model] = [grid] 
        try:
            key = met_dict[model]
            if met not in key:                                                                                                                                                                                                           
                key.append(met) 
                met_dict[model] = key
        except:
            met_dict[model] = [met]     

    try:
        key = model_dict[model]
        if year_range not in key:
            key.append(year_range) 
            model_dict[model] = key
    except:
        model_dict[model] = [year_range]   
    
#get model and date range
model_name = raw_input('\nChoose Model.\n%s\n'%('   '.join(i for i in model_dict)))
model_range = model_dict[model_name]
if len(model_range) == 1:
    model_range = model_range[0]
else:
    model_range = raw_input('\nChoose Date Range.\n%s\n'%('   '.join(i for i in model_range)))

#get model version
valid_versions = []
try:
    versions = version_dict[model_name]
    for mf in model_files:
        for i in versions:
            test_string = '%s_SURFACE_%s_%s_%s'%(model_name,species,model_range,i)
            if test_string in mf:
                valid_versions.append(i)     
    valid_versions = set(valid_versions)
    valid_versions = [i for i in valid_versions]   
except:
    valid_versions = ''
if valid_versions != '':
    if len(valid_versions) > 1:
        version = raw_input('\nChoose Model Version.\n%s\n'%('   '.join(i for i in valid_versions)))
    else:
        version = valid_versions[0]

#get grid version
valid_grids = []
try:
    grids = grid_dict[model_name]
    for mf in model_files:
        for i in grids:
            test_string = '%s_SURFACE_%s_%s_%s_%s'%(model_name,species,model_range,version,i)
            if test_string in mf:
                valid_grids.append(i)
    valid_grids = set(valid_grids)
    valid_grids = [i for i in valid_grids]
except:
    valid_grids = ''
if valid_grids != '':
    if len(valid_grids) > 1:
        grid = raw_input('\nChoose Model Grid.\n%s\n'%('   '.join(i for i in valid_grids)))
    else:
        grid = valid_grids[0]

#get met version
valid_mets = []
try:
    mets = met_dict[model_name]
    for mf in model_files:
        for i in mets:
            test_string = '%s_SURFACE_%s_%s_%s_%s_%s'%(model_name,species,model_range,version,grid,i)
            if test_string in mf:
                valid_mets.append(i)
    valid_mets = set(valid_mets)
    valid_mets = [i for i in valid_mets]
except:
    valid_mets = ''
if valid_mets != '':
    if len(valid_mets) > 1:
        met = raw_input('\nChoose Model Meteorology.\n%s\n'%('   '.join(i for i in valid_mets)))
    else:
        met = valid_mets[0]

#put together model file name
if (valid_versions == '') & (valid_grids == '') & (valid_mets == ''):
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'.nc'
    version,grid,met = '','',''
elif (valid_versions != '') & (valid_grids == '') & (valid_mets == ''):
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'.nc'
    grid,met = '',''
elif (valid_versions != '') & (valid_grids != '') & (valid_mets == ''):
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'_'+grid+'.nc'
    met = ''
elif (valid_versions != '') & (valid_grids != '') & (valid_mets != ''):
    model_f = '/work/home/db876/plotting_tools/model_files/'+model_name+'_'+'SURFACE'+'_'+species+'_'+model_range+'_'+version+'_'+grid+'_'+met+'.nc'

root_grp = Dataset(model_f)
model_var = root_grp.variables[species.lower()][:]
model_date = root_grp.variables['date'][:]
model_time = root_grp.variables['time'][:]
lat_e = root_grp.variables['lat_edges'][:]
lon_e = root_grp.variables['lon_edges'][:]
lat_c = root_grp.variables['lat_centre'][:]
lon_c = root_grp.variables['lon_centre'][:]
grid_size = root_grp.variables['grid_size'][:]
grid_size = grid_size[0]

gridbox_count = len(lat_c)*len(lon_c)

#process obs
#----------------------------------------
#data range

#read in obs data
try:
    root_grp = Dataset('../process/GLOBAL_SURFACE_%s_%s_%s.nc'%(species,model_range[:4],model_range[-4:]))
    note = ''
except:
    all_obs_files = glob.glob('../process/GLOBAL_SURFACE_%s*'%(species)) 
    valid_obs_files = []
    for f in all_obs_files:
        root_grp = Dataset(f)
        valid_refs_dict = root_grp.groups
        valid_refs = []
        for i in valid_refs_dict.keys():
            valid_refs = np.append(valid_refs,i)
        valid_ref = valid_refs[0]
        site_group = root_grp.groups[valid_ref]
        test_obs_var = site_group.variables[species.lower()][:]
        if len(model_var) == len(test_obs_var):
            valid_obs_files.append(f)
    print valid_obs_files
    
    if len(valid_obs_files) == 1:
        m_start_year, m_end_year = start_year,end_year
        file = valid_obs_files[0]
        start_year = file[-12:-8]
        end_year = file[-7:-3]
        note = '- BE AWARE OBS. IS %s to %s,  MODEL IS %s to %s'%(start_year,end_year,m_start_year,m_end_year) 
        print '\nNo equivalent observational data range. Using range between %s and %s that matches in shape.\n'%(start_year,end_year)    
        root_grp = Dataset(file)        

    elif len(valid_obs_files) > 1:
        print '\nNo equivalent observational data range. Choose a range that matches in shape to compare against.\n' 
        s_years = []
        e_years = []
        for i in range(len(valid_files)):
            file = valid_files[i]
            s_years.append(file[-12:-8])
            e_years.append(file[-7:-3])
        ranges = ['%s_%s'%(a,b) for a,b in zip(s_years,e_years)] 
        chosen_range = raw_input('\nChoose Range.\n%s\n'%('   '.join(i for i in ranges)))
        chosen_file = valid_files[ranges.index(chosen_range)]
        m_start_year, m_end_year = start_year,end_year
        start_year = int(chosen_range[:4])
        end_year = int(chosen_range[5:])
        note = '- BE AWARE OBS. IS %s to %s,  MODEL IS %s to %s'%(start_year,end_year,m_start_year,m_end_year)
        root_grp = Dataset(chosen_file)
    else:
        print '\nNo equivalent observational data range. And no observational data range that matches in shape. Adios.\n'
        1+'a'

valid_refs_dict = root_grp.groups
valid_refs = []
for i in valid_refs_dict.keys():
    valid_refs = np.append(valid_refs,i)

site_ref = raw_input('\nChoose site from list. Sites with full set of yearly files between %s & %s are:\n%s\n'%(start_year,end_year,'   '.join(i for i in valid_refs)))

#read in specific site data
site_group = root_grp.groups[site_ref]

#read in variables for site
obs_var = site_group.variables[species.lower()][:]
obs_date = site_group.variables['date'][:]
obs_time = site_group.variables['time'][:]
obs_lat = site_group.latitude
obs_lon = site_group.longitude
obs_alt = site_group.altitude
obs_group = site_group.process_group

obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
obs_ave = np.ma.average(obs_var_mask)

#if using a different obs time rnage to compare model against, take the obs date and time to be the model's also.
if note != '':
    model_date = obs_date[:]
    model_time = obs_time[:]

#----------------------------------------
#find model data gridbox to compare with obs.

#get model gridbox for obs site
lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)

model_var = model_var[:,lat_n,lon_n]
model_var = model_var*1e9

model_var_mask = np.ma.masked_where(model_var<=0,model_var)
model_ave = np.ma.average(model_var_mask)

#--------------------------------------------
#take half daily average of obs and model

obs_time = modules.date_process(obs_date,obs_time,start_year)
model_time = modules.date_process(model_date,model_time,start_year)

divisor = 6

#take half daily average of obs
total_len = len(obs_var_mask)/divisor
start = 0
end = divisor
ave_obs_var = []
ave_obs_time = []
for i in range(total_len):
    ave = np.ma.average(obs_var_mask[start:end])
    ave_obs_time=np.append(ave_obs_time,obs_time[start])
    ave_obs_var=np.append(ave_obs_var,ave)
    start+=divisor
    end+=divisor
 
obs_test = ave_obs_var >= 0
ave_obs_var = ave_obs_var[obs_test]
ave_obs_time = ave_obs_time[obs_test]

print ave_obs_time

#take half daily average of model
total_len = len(model_var_mask)/divisor
start = 0
end = divisor
ave_model_var = []
ave_model_time = []
for i in range(total_len):
    ave = np.ma.average(model_var_mask[start:end])
    ave_model_time=np.append(ave_model_time,model_time[start])
    ave_model_var=np.append(ave_model_var,ave)
    start+=divisor
    end+=divisor

model_test = ave_model_var >= 0
ave_model_var = ave_model_var[model_test]
ave_model_time = ave_model_time[model_test]

#take lomb for obs and model
mctest = True
nsim = int(raw_input('\nChoose Number of Redfit Calculation iterations.\n'))
n50 = int(raw_input('\nChoose number of WOSA segments. (1 for standard LSP)\n'))
ofac = int(raw_input('\nChoose oversampling factor. (Typically 4.)\n'))

#obs_periods,obs_mag,obs_fr,obs,fi,obs_red_periods,obs_red_mag,obs_gredth,obs_fac80,obs_fac85,obs_fac90,obs_fac95,obs_fac99,obs_fac99_9,obs_faccrit,obs_tau,obs_corr = redfit.redcore(nsim,n50,mctest,ave_obs_time,ave_obs_var)
model_periods,model_mag,model_fr,model_fi,model_red_periods,model_red_mag,model_gredth,model_fac80,model_fac85,model_fac90,model_fac95,model_fac99,model_fac99_9,model_faccrit,model_tau,model_corr = redfit.redcore(nsim,n50,mctest,ave_model_time,ave_model_var)

print 'model fr len = ',  len(model_fr)
#-------------------------------------------------------------------------------
#Do IFFT of altered spectra -removing significant periods, and one unaltered.
#altered spectra provides red noise est baseline

#IFFT for unaltered spectra

##use ifft to get time series back from unadjusted spectra
#complex Fourier spectrum which corresponds to the Lomb-Scargle periodogram:
F = [0]*((len(model_fr)*2)+1)
#set first real value to average 
F[0] = complex(model_ave,0)

#Get reverse real and imaginary values
rev_model_fr=np.copy(model_fr[::-1])
rev_model_fi=np.copy(model_fi[::-1])

#Fill Fourier Spectrum real and imaginary values
for i in range(len(model_fr)):
    F[i+1] = complex(model_fr[i],model_fi[i])
for i in range(len(model_fr),len(model_fr)*2):
    F[i+1] = complex(rev_model_fr[i-len(model_fr)],-rev_model_fi[i-len(model_fr)])
F = np.array(F)    
    
print 'len F = ', len(F)    

#Take ifft and just take real values
model_ts = numpy.fft.ifft(F)
model_ts = model_ts.astype('float64')

print 'ts len = ', len(model_ts)

time = np.arange(0,len(model_ts)*1.,0.25)

actual_model_ts_len = (len(model_ts)/ofac) + np.mod(len(model_ts),ofac)

time = time[-actual_model_ts_len:]
model_ts = model_ts[-actual_model_ts_len:]

#model_ts=model_ts/scipy.signal.hanning(len(model_ts))

time = time-time[0]

print time

#remove first 1/6th and last 1/6th because of damage due to windowing
#wind_cut = len(model_ts)/6
#time = time[wind_cut:-wind_cut]
#model_ts = model_ts[wind_cut:-wind_cut]

model_ts = model_ts+model_ave

print 'ts len = ', len(model_ts)

print np.min(model_ts)
#get spectra from unaltered ts
model_periods_2,model_mag_2,model_fr_2,model_fi_2,red_periods_2,red_mag_2,gredth_2,fac80_2,fac85_2,fac90_2,fac95_2,fac99_2,fac99_9_2,faccrit_2,tau_2,corr_2 = redfit.redcore(nsim,n50,mctest,time,model_ts)

#IFFT for altered spectra
#Take out peaks > 95 % significance

model_sig95 = model_gredth*model_fac95
model_sig_test = np.array(model_mag >= model_sig95) 
model_sig_inds = np.where(model_sig_test == True)
model_sig_inds = [i for i in model_sig_inds[0]]

#Where values are >= significance level, change real and imaginary numbers to interpolate between surrounding points.
#also remove sidelobes of peaks to certain % of peak max.

print 'model fr len = ', len(model_fr)

alter_model_mag,alter_model_fr,alter_model_fi,sig_sidelobe_inds = redfit.sidelobe_peak_remove(model_mag,model_fr,model_fi,model_sig_inds,50,model_periods)

print 'model fr len = ', len(alter_model_fr)

##use ifft to get time series back from adjusted spectra
#complex Fourier spectrum which corresponds to the Lomb-Scargle periodogram:
F = [0]*((len(alter_model_fr)*2)+1)
#set first real value to average 
F[0] = complex(model_ave,0)

#Get reverse real and imaginary values
rev_alter_model_fr=np.copy(alter_model_fr[::-1])
rev_alter_model_fi=np.copy(alter_model_fi[::-1])

#Fill Fourier Spectrum real and imaginary values
for i in range(len(alter_model_fr)):
    F[i+1] = complex(alter_model_fr[i],alter_model_fi[i])
for i in range(len(alter_model_fr),len(alter_model_fr)*2):
    F[i+1] = complex(rev_alter_model_fr[i-len(alter_model_fr)],-rev_alter_model_fi[i-len(alter_model_fr)])
F = np.array(F)    

print 'len F = ', len(F)

#Take ifft and just take real values
alter_model_ts = numpy.fft.ifft(F)
alter_model_ts = alter_model_ts.astype('float64')

print 'ts len = ', len(alter_model_ts)

print alter_model_ts

alter_time = np.arange(0,len(alter_model_ts)*1.,0.25)

actual_model_ts_len = (len(alter_model_ts)/ofac) + np.mod(len(alter_model_ts),ofac)

alter_time = alter_time[-actual_model_ts_len:]
alter_model_ts = alter_model_ts[-actual_model_ts_len:]

#alter_model_ts=alter_model_ts/scipy.signal.hanning(len(alter_model_ts))

alter_time = alter_time-alter_time[0]

#remove first 1/6th and last 1/6th because of damage due to windowing
#wind_cut = len(alter_model_ts)/6
#alter_time = alter_time[wind_cut:-wind_cut]
#alter_model_ts = alter_model_ts[wind_cut:-wind_cut]

alter_model_ts = alter_model_ts+model_ave

print 'ts len = ', len(alter_model_ts)

#get rednoise estimation from altered ts
print np.min(alter_model_ts)
print alter_model_ts
print len(alter_model_ts)
model_periods_3,model_mag_3,model_fr_3,model_fi_3,red_periods_3,red_mag_3,gredth_3,fac80_3,fac85_3,fac90_3,fac95_3,fac99_3,fac99_9_3,faccrit_3,tau_3,corr_3 = redfit.redcore(nsim,n50,mctest,alter_time,alter_model_ts)

plt.plot(ave_model_var)
plt.plot(model_ts)
plt.plot(alter_model_ts)
plt.show()

#plt.plot(ave_model_time,ave_model_var,label='orig ts')
#plt.plot(time,model_ts,label = 'ts')


# #put old mags back in where significant peaks were removed, not sidelobes. As length of spectra has changed need to alter closest point to where original period was.
# #if indices are together, then points are same as same peak, so just get peak with biggest mag
# 
# alter_x_mag = np.copy(x_mag)
# 
# same_peak = []
# key_peak_inds = []
# for i in range(len(sig_sidelobe_inds)):
#     try:
#         diff = sig_sidelobe_inds[i+1] - sig_sidelobe_inds[i]
#     except:
#         diff = 0
#     if diff == 1:
#         same_peak=np.append(same_peak,sig_sidelobe_inds[i])
#     else:
#         same_peak=np.append(same_peak,sig_sidelobe_inds[i])
#         if len(same_peak) > 1:
#             same_peak = same_peak.astype(int)
#             max_ind = np.argmax(original_model_mag[same_peak])
#             max_ind = same_peak[max_ind]
#             key_peak_inds = np.append(key_peak_inds,max_ind)
#         else:
#             same_peak = same_peak.astype(int)
#             key_peak_inds = np.append(key_peak_inds,sig_sidelobe_inds[i])
#             
#         same_peak = []
#     
# key_peak_inds = key_peak_inds.astype(int)
# original_peak_periods = model_periods[key_peak_inds]
# original_peak_mags  = original_model_mag[key_peak_inds]
# 
# for i in range(len(original_peak_periods)):
#     peak_period = original_peak_periods[i]
#     peak_mag = original_peak_mags[i]
#     new_peak_ind = min(range(len(x_periods)), key=lambda i: abs(x_periods[i]-peak_period))
#     x_mag[new_peak_ind] = peak_mag
# 
# sig95 = new_gredth*new_fac95
# sig99 = new_gredth*new_fac99
# sig_inds = []
# sig_periods = []
# for i in range(len(new_gredth)):
#     if x_mag[i] >= (new_gredth[i]*new_fac99):
#         sig_inds.append(i)
#         sig_periods.append(x_periods[i])
# 
# print 'Sig periods = ', sig_periods



# #if sig inds are next to each other, just leave 1 peak in array by deleting other points
# same_peak = []
# 
# for i in range(len(sig_inds)):
#     try:
#         diff = sig_inds[i+1] - sig_inds[i]
#     except:
#         diff = 0
#     if diff == 1:
#         same_peak=np.append(same_peak,sig_inds[i])
#     else:
#         same_peak=np.append(same_peak,sig_inds[i])
#         if len(same_peak) > 1:
#             same_peak = same_peak.astype(int)
#             max_ind = np.argmax(original_model_mag[same_peak])
#             max_ind = same_peak[max_ind]
#             
#             not_max_inds = []
#             for j in same_peak:
#                 if j != max_ind:
#                     not_max_inds.append(j)
#             
#             x_mag = np.delete(x_mag,not_max_inds)
#             x_periods = np.delete(x_periods,not_max_inds)
#              
#         same_peak = []
 
#----------------------------------------------
#plotting

#set plotting area & background to white
fig=plt.figure(figsize=(23,12.3))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1)

#plt.loglog(obs_periods,obs_mag, color='black', label = 'Obs.')
#plt.loglog(obs_red_periods,obs_gredth,color='yellow')
#plt.loglog(obs_red_periods,obs_gredth*obs_fac95,color='pink', label=' Obs. 95% Significance')
#plt.loglog(obs_red_periods,obs_gredth*obs_fac99,color='orange', label=' Obs. 99% Significance')


plt.loglog(model_periods_2,model_mag_2, marker='x', color='red', label = '%s'%(model_name))
#plt.loglog(model_periods_3,model_mag_3, marker='x', color='red', label = '%s'%(model_name))
#plt.loglog(model_periods,model_mag, marker='x', color='red', label = '%s'%(model_name))

plt.semilogx(red_periods_3,gredth_3,color='blue')
plt.semilogx(red_periods_3,gredth_3*fac95_3,color='green', label='Model 95% Significance')
plt.semilogx(red_periods_3,gredth_3*fac99_3,color='purple', label='Model 99% Significance')
plt.semilogx(red_periods_3,gredth_3*fac99_9_3,color='black', label='Model 99.9% Significance')

#plt.loglog(model_red_periods,model_gredth,color='blue')
#plt.loglog(model_red_periods,model_gredth*model_fac95,color='green', label='Model 95% Significance')
#plt.loglog(model_red_periods,model_gredth*model_fac99,color='purple', label='Model 99% Significance')
        
def form2(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.2f' % x	

def form5(x, pos):
	""" This function returns a string with 3 decimal places, given the input x"""
	return '%.5f' % x

xformatter = FuncFormatter(form2)
yformatter = FuncFormatter(form5)

plt.grid(True)
leg=plt.legend(loc=4, prop={'size':24})
leg.get_frame().set_alpha(0.4)
ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))
plt.xlabel('Period (Days)',fontsize = 21)
plt.ylabel('Amplitude (ppb)',fontsize = 21)

ax.tick_params(axis='both', which='major', labelsize=18)
for tick in ax.get_xaxis().get_major_ticks():
    tick.set_pad(11.)
fig.tight_layout()
plt.show()
