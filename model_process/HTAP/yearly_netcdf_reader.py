from scipy.io import netcdf
import numpy as np
import glob
import matplotlib.pyplot as plt
import lomb
import nappy

#Set desired location
location = 'Cape Verde'

#List of locations

#Cape Point
#Cape Verde
#Mace Head
#Mauna Loa
#Ragged Point
#Tudor Hill
#Tutuila


#Set desired species
species = 'NO2'

#List of species

#O3
#NO2

def movingaverage(interval, window_size):
    window= np.ones(int(window_size))/float(window_size)
    
    return np.convolve(interval, window, 'same')

def read_model(files):
    variables =[]
    time = []

    for i in files:
        read = netcdf.netcdf_file(i, 'r')
#Extract index of latitude for relevant location
        mod_lat = abs(read.variables['lat'][:]) - real_lat
        f = lambda lst, target: lst[np.argmin(np.abs(np.array(lst) - target))]
        closest_lat = f(mod_lat,0)

        list_mod_lat = np.array(mod_lat).tolist()
        round_mod_lat = []
    #for i in list_mod_lat:
    #round_mod_lat.append(round(i,3))
        lat_index = list_mod_lat.index(closest_lat)

    #Extract index of longitude for relevant location
        mod_lon = abs(read.variables['lon'][:]) - real_lon
        f = lambda lst, target: lst[np.argmin(np.abs(np.array(lst) - target))]
        closest_lon = f(mod_lon,0)

        list_mod_lon = np.array(mod_lon).tolist()
        round_mod_lon = []
    #for i in list_mod_lon:
    #round_mod_lon.append(round(i,3))
        lon_index = list_mod_lon.index(closest_lon)
    
        variable_array = read.variables[spec][:,lon_index,lat_index]
        time_array = read.variables['time'][:]
        time.append(time_array)
        variables.append(variable_array)

    variables = [item for sublist in variables for item in sublist]
    time = [item for sublist in time for item in sublist]

    variables = np.array(variables)
    variables = variables[:]*1e12
    standard_deviation_model_p = np.std(variables)
    mean_model_p = np.mean(variables)
    normal_model = variables-mean_model_p
    normal_model = normal_model/standard_deviation_model_p

    time =  np.array(time)

    return time, variables, normal_model

#Input Latitude and Longitude for relevant location
if location == 'Cape Point':
    real_lat = -34.35
    real_lon = 18.48

if location == 'Cape Verde':
    real_lat = 16.848   #For CVO
    real_lon = -24.871

if location == 'Mace Head':
    real_lat = 53.33
    real_lon = -9.90

if location == 'Mauna Loa':
    real_lat = 19.54
    real_lon = -155.58

if location == 'Ragged Point':
    real_lat = 13.17
    real_lon = -59.43

if location == 'Tudor Hill':
    real_lat = 32.27
    real_lon = -64.87

if location == 'Tutuila':
    real_lat = -14.24
    real_lon = -170.57

if species == 'O3':
    spec = 'vmr_o3'

if species == 'NO2':
    spec = 'vmr_no2'

#Read model files into memory
camchem_3311_time, camchem_3311_o3, norm_camchem_3311_o3 =  read_model(glob.glob('model_files/CAMCHEM-3311m13*'))
camchem_3514_time, camchem_3514_o3, norm_camchem_3514_o3 =  read_model(glob.glob('model_files/CAMCHEM-3514*'))
cam_sr1_time, cam_sr1_o3, norm_cam_sr1_o3 =  read_model(glob.glob('model_files/CAM_SR1*'))
chaser_time, chaser_o3, norm_chaser_o3 =  read_model(glob.glob('model_files/CHASER*'))
frsgcuci_time, frsgcuci_o3, norm_frsgcuci_o3 =  read_model(glob.glob('model_files/FRSGCUCI*'))
gemaq_time, gemaq_o3, norm_gemaq_o3 =  read_model(glob.glob('model_files/GEMAQ*'))
geoschem_time, geoschem_o3, norm_geoschem_o3 =  read_model(glob.glob('model_files/GEOSChem*'))
giss_time, giss_o3, norm_giss_o3 =  read_model(glob.glob('model_files/GISS-PUCCINI-modelEaer_SR1_sfc*'))
giss_alt_time, giss_alt_o3, norm_giss_alt_o3 =  read_model(glob.glob('model_files/GISS_PUCCINI_modelE_alt_SR1*'))
inca_time, inca_o3, norm_inca_o3 =  read_model(glob.glob('model_files/INCA*'))
llnl_time, llnl_o3, norm_llnl_o3 =  read_model(glob.glob('model_files/LLNL*'))
mozart_time, mozart_o3, norm_mozart_o3 =  read_model(glob.glob('model_files/MOZARTGFDL*'))
mozech_time, mozech_o3, norm_mozech_o3 =  read_model(glob.glob('model_files/MOZECH*'))
oslo_time, oslo_o3, norm_oslo_o3 =  read_model(glob.glob('model_files/OsloCTM2*'))
#tm5_time, tm5_o3, norm_tm5_o3 =  read_model(glob.glob('model_files/TM5-JRC*'))


#Read in obs
#now read in the observations

myfile=nappy.openNAFile('York_merge_Cape_verde_1hr_R1.na')
myfile.readData()

#ppy.openNAFile('York_merge_Cape_verde_1hr_R1.na')
k_var1=myfile["VNAME"].index('Ozone mixing ratio (ppbV)_(Mean)')

# OK need to conver values from a list to a numpy array

time=np.array(myfile['X'])
var1=np.array(myfile['V'][k_var1])
valids1=var1 > 0
time=time[valids1]
var=var1[valids1]

valids2= time <= 730
obs_time=time[valids2]
var2=var[valids2]

valids3= obs_time >= 365
obs_time=obs_time[valids3]
var2=var2[valids3]

standard_deviation_obs_p = np.std(var2)
mean_obs_p = np.mean(var2)
obs_var = var2-mean_obs_p
obs_var = obs_var/standard_deviation_obs_p

standard_deviation_obs_p = np.std(var)
mean_obs_p = np.mean(var)
obs_var_full = var-mean_obs_p
obs_var_full = obs_var_full/standard_deviation_obs_p

#Plot them all up. 
fig=plt.figure(figsize=(20,12))

#Plot up standard conc. v time plot
ax1= fig.add_subplot(2, 1, 1)
fig.subplots_adjust(hspace=0.3)
#plt.plot(obs_time, var2 , color='k', label='2007 CV Obs. ')


#plt.plot(obs_time, var2, color='k', label='2007 CV Obs. ')
plt.plot(camchem_3311_time, camchem_3311_o3, color='b', label='CAMCHEM-3311m13')
plt.plot(camchem_3514_time, camchem_3514_o3, color='darkgreen', label='CAMCHEM-3514')
plt.plot(cam_sr1_time, cam_sr1_o3, color='c', label='CAM_SR1')
plt.plot(chaser_time, chaser_o3, color='g', label='CHASER')
plt.plot(frsgcuci_time, frsgcuci_o3, color='lightslategrey', label='FRSGCUCI')
plt.plot(gemaq_time, gemaq_o3, color='m', label='GEMAQ')
plt.plot(geoschem_time, geoschem_o3, color='r', label='GEOSChem')
plt.plot(giss_time, giss_o3, color='y', label='GISS')
plt.plot(giss_alt_time, giss_alt_o3, color='orangered', label='GISS_alt')
plt.plot(inca_time,inca_o3, color='saddlebrown', label ='INCA')
plt.plot(llnl_time, llnl_o3, color='lightskyblue', label ='LLNL')
plt.plot(mozart_time,mozart_o3, color='gold', label = 'MOZARTGFDL')
plt.plot(mozech_time,mozech_o3, color='tomato', label ='MOZECH')
plt.plot(oslo_time,oslo_o3, color='greenyellow', label ='OSLO')
plt.grid(True)
#plt.xlim(300,301)
#plt.ylim(0,250)
leg=plt.legend(loc=1)
leg.get_frame().set_alpha(0.4)
plt.xlabel('Time (Days since 2001)')
plt.ylabel('%s (%s)' % ('Conc.','ppbv'))
plt.title('HTAP Model %s V Time at %s'% (species,location))
 
#Define sampling frequency
samp_freq = 24

#Lomb-scargle plot
ax3= fig.add_subplot(2, 1, 2)

#Plot axis period lines and labels
annotate_line_y=np.arange(1e-10,1e4,1)
freq_year = [345]*len(annotate_line_y)
plt.plot(freq_year, annotate_line_y,'r--',alpha=0.4)
plt.text(345, 1e-7, '1 Year', fontweight='bold')

#Model lomb
freq_obs, power_obs, nout, jmax, prob = lomb.fasper(obs_time, obs_var, 1., samp_freq)
freq_obs_full, power_obs_full, nout, jmax, prob = lomb.fasper(time, obs_var_full, 1., samp_freq)
freq_camchem3311, power_camchem3311, nout, jmax, prob2 = lomb.fasper(camchem_3311_time,norm_camchem_3311_o3, 1., samp_freq)
freq_camchem3514, power_camchem3514, nout, jmax, prob2 = lomb.fasper(camchem_3514_time,norm_camchem_3514_o3, 1., samp_freq)
freq_cam_sr1, power_cam_sr1, nout, jmax, prob2 = lomb.fasper(cam_sr1_time,norm_cam_sr1_o3, 1., samp_freq)
freq_chaser, power_chaser, nout, jmax, prob2 = lomb.fasper(chaser_time,norm_chaser_o3, 1., samp_freq)
freq_frsgcuci, power_frsgcuci, nout, jmax, prob2 = lomb.fasper(frsgcuci_time,norm_frsgcuci_o3, 1., samp_freq)
freq_gemaq, power_gemaq, nout, jmax, prob2 = lomb.fasper(gemaq_time,norm_gemaq_o3, 1., samp_freq)
freq_geoschem, power_geoschem, nout, jmax, prob2 = lomb.fasper(geoschem_time,norm_geoschem_o3, 1., samp_freq)
freq_giss, power_giss, nout, jmax, prob2 = lomb.fasper(giss_time,norm_giss_o3, 1., samp_freq)
freq_giss_alt, power_giss_alt, nout, jmax, prob2 = lomb.fasper(giss_alt_time,norm_giss_alt_o3, 1., samp_freq)
freq_inca, power_inca, nout, jmax, prob2 = lomb.fasper(inca_time,norm_inca_o3, 1., samp_freq)
freq_llnl, power_llnl, nout, jmax, prob2 = lomb.fasper(llnl_time,norm_llnl_o3, 2.0001, samp_freq)
freq_mozart, power_mozart, nout, jmax, prob2 = lomb.fasper(mozart_time,norm_mozart_o3, 1., samp_freq)
freq_mozech, power_mozech, nout, jmax, prob2 = lomb.fasper(mozech_time,norm_mozech_o3, 1., samp_freq)
freq_oslo, power_oslo, nout, jmax, prob2 = lomb.fasper(oslo_time,norm_oslo_o3, 1., samp_freq)
#freq_tm5, power_tm5, nout, jmax, prob2 = lomb_edit.fasper(tm5_time,norm_tm5_o3, 1., samp_freq)

#Divide output by sampling 
power_obs = power_obs/samp_freq
power_obs_full = power_obs_full/samp_freq
power_camchem3311 = power_camchem3311/samp_freq
power_camchem3514 = power_camchem3514/samp_freq
power_cam_sr1 = power_cam_sr1/samp_freq
power_chaser = power_chaser/samp_freq
power_frsgcuci = power_frsgcuci/samp_freq
power_gemaq = power_gemaq/samp_freq
power_geoschem = power_geoschem/samp_freq
power_giss = power_giss/samp_freq
power_giss_alt = power_giss_alt/samp_freq
power_inca = power_inca/samp_freq
power_llnl = power_llnl/samp_freq
power_mozart = power_mozart/samp_freq
power_mozech = power_mozech/samp_freq
power_oslo = power_oslo/samp_freq
#power_tm5 = power_tm5/samp_freq


#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
    #nyquist_freq_lomb_model = frequencies[-1]
    #Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
    #print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 

#Gaussian Smoothing, set window averaging size
window_size = 4

power_obs = movingaverage(power_obs,window_size)
power_obs_full = movingaverage(power_obs_full,window_size)
power_camchem3311 = movingaverage(power_camchem3311,window_size)
power_camchem3514 = movingaverage(power_camchem3514,window_size)
power_cam_sr1 = movingaverage(power_cam_sr1,window_size)
power_chaser = movingaverage(power_chaser,window_size)
power_frsgcuci = movingaverage(power_frsgcuci,window_size)
power_gemaq = movingaverage(power_gemaq,window_size)
power_geoschem = movingaverage(power_geoschem,window_size)
power_giss = movingaverage(power_giss,window_size)
power_giss_alt = movingaverage(power_giss_alt,window_size)
power_inca = movingaverage(power_inca,window_size)
power_llnl = movingaverage(power_llnl,window_size)
power_mozart = movingaverage(power_mozart,window_size)
power_mozech = movingaverage(power_mozech,window_size)
power_oslo = movingaverage(power_oslo,window_size)
#power_tm5 = movingaverage(power_tm5,10)

#plot up
plt.loglog(1./freq_obs, power_obs, 'kx', linestyle ='-', alpha = 0.75,markersize=2, label='2007 CV Obs. ')
plt.loglog(1./freq_obs_full, power_obs_full, color='teal', marker='x', linestyle ='-', alpha = 0.75,markersize=2, label='2006-2012 CV Obs. ')
plt.loglog(1./freq_camchem3311, power_camchem3311, 'bx', linestyle ='-', alpha = 0.75,markersize=2, label='CAMCHEM-3311m1 ')
plt.loglog(1./freq_camchem3514, power_camchem3514, color='darkgreen', marker='x', linestyle ='-', alpha = 0.75,markersize=2, label='CAMCHEM-3514 ')
plt.loglog(1./freq_cam_sr1, power_cam_sr1, 'cx', linestyle ='-', alpha = 0.75,markersize=2, label='CAM_SR1 ')
plt.loglog(1./freq_chaser, power_chaser, 'gx', linestyle ='-', alpha = 0.75,markersize=2, label='CHASER ')
plt.loglog(1./freq_frsgcuci, power_frsgcuci, color='lightslategrey', marker ='x', linestyle ='-', alpha = 0.75,markersize=2, label='FRSGCUCI ')
plt.loglog(1./freq_gemaq, power_gemaq, 'mx', linestyle ='-', alpha = 0.75,markersize=2, label='GEMAQ ')
plt.loglog(1./freq_geoschem, power_geoschem, 'rx', linestyle ='-', alpha = 0.75,markersize=2, label='GEOSChem ')
plt.loglog(1./freq_giss, power_giss, 'yx', linestyle ='-', alpha = 0.75,markersize=2, label='GISS ')
plt.loglog(1./freq_giss_alt, power_giss_alt, color='orangered', marker ='x',linestyle ='-', alpha = 0.75,markersize=2, label='GISS_alt ')
plt.loglog(1./freq_inca, power_inca, color='saddlebrown', marker ='x',linestyle ='-', alpha = 0.75,markersize=2, label='INCA ')
plt.loglog(1./freq_llnl, power_llnl, color='lightskyblue', marker ='x',linestyle ='-', alpha = 0.75,markersize=2, label='LLNL ')
plt.loglog(1./freq_mozart, power_mozart, color='gold', marker ='x',linestyle ='-', alpha = 0.75,markersize=2, label='MOZART ')
plt.loglog(1./freq_mozech, power_mozech, color='tomato', marker ='x',linestyle ='-', alpha = 0.75,markersize=2, label='MOZECH ')
plt.loglog(1./freq_oslo, power_oslo, color='greenyellow', marker ='x',linestyle ='-', alpha = 0.75,markersize=2, label='OSLO ')
#plt.loglog(1./freq_tm5, power_tm5, color='orangered', marker ='x',linestyle ='-', alpha = 0.75,markersize=2, label='TM5-JRC ')
    #make index for high frequency measure of obs. and model PSD. Say is between min period(2 hours) and 1 day
#model_periods = 1./fx

#model_periods = [x for x in model_periods if x <= 1]

    #cut psds to same length as period   
#psd_model=fy[len(model_periods):]

#ave_mod_per = np.average(psd_model)

plt.grid(True)
leg=plt.legend(loc=7)
#leg.get_frame().set_alpha(0.4)
#plt.text(1e-2, 1e2,'Model ave PSD from 2hrs:1day = %4.10f' %(ave_mod_per),  fontweight='bold')
#plt.xlim(343, 346)
plt.ylim(1e-1,1e3)
plt.xlabel('Period (Days)')
plt.ylabel(r'PSD $(ppb^{2}/days^{-1})$')
plt.title('Lomb-Scargle %s Power V Period at %s' % (species,location))

plt.show()
