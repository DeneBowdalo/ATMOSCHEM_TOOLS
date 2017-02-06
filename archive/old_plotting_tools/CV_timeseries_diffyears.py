import numpy as np
import csv
import glob
import datetime as datetime 
import matplotlib.pyplot as plt
import logging as log
import nappy
import scipy.stats as stats
import lomb
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq, rfft

#List of plotable species:

#O3
#CO
#NO
#NO2
#C2H6
#C3H8
#DMS
#ISOPRENE
#ACETONE
#TEMP
#SURFACE_PRES
#WINDSPEED
#SURFACE_SOLAR_RADIATION
#ABS_HUMIDITY
#REL_HUMIDITY

#Input species below from list above

species = 'O3'


#Gives species exact model tags for convenience
if species == 'ISOPRENE':
    species = 'TRA_6'

elif species == 'ACETONE':
    species = 'ACET'

elif species == 'TEMP':
    species = 'GMAO_TEMP'

elif species == 'SURFACE_PRES':
    species = 'GMAO_PSFC'

elif species == 'WINDSPEED':
    species = 'GMAO_WIND'

elif species == 'SURFACE_SOLAR_RADIATION':
    species = 'GMAO_RADSW'

elif species == 'ABS_HUMIDITY':
    species = 'GMAO_ABSH'

elif species == 'REL_HUMIDITY':
    species = 'GMAO_RHUM'


#make index for high frequency measure of obs. and model PSD. Say is between min period(2 hours) and 1 day 
def period_percent_diff(start, end, fb, fy, obs_periods, model_periods):

    cut_obs_periods = [x for x in obs_periods if start <= x < end]
    cut_model_periods = [x for x in model_periods if start <= x < end]

    #cut psds to same length as period   
    psd_obs = fb[len(cut_obs_periods):]
    psd_model=fy[len(cut_model_periods):]

    ave_obs_per = np.average(psd_obs)
    ave_mod_per = np.average(psd_model)

    if ave_obs_per > ave_mod_per:
        biggest = 'Obs.'
        power_diff = ave_obs_per - ave_mod_per
        average = (ave_mod_per + ave_obs_per) / 2
        percent = (power_diff / average) * 100
    elif ave_mod_per > ave_obs_per:
        biggest = 'Model'
        power_diff = ave_mod_per - ave_obs_per
        average = (ave_mod_per + ave_obs_per) / 2
        percent = (power_diff / average) * 100
        
    return percent 

                                     
#reads in species name and processes unique variables
def plot(species):

    #Set model_cut switch to default 0, if want to do more complicated cuts from model field, specify model_cut_switch == 1 in species definitions
     #Vice versa with obs_switch
    model_cut_switch = 0
    obs_switch = 0
    ofac = 1
    if species == 'O3':
        units = 'ppbV'
        first_label_pos = 3
        obs_data_name = 'Ozone mixing ratio (ppbV)_(Mean)'
        unit_cut= 1e9
        species_type = 'Conc.'
        actual_species_name = 'O3'

    elif species == 'CO':
        units = 'ppbV'
        first_label_pos = 1
        obs_data_name = 'CO mixing ratio (ppbV)_(Mean)'         
        unit_cut= 1e9
        species_type = 'Conc.'
        actual_species_name = 'CO'
        ofac = 2.0001

    elif species == 'NO':
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'NO mixing ratio (pptv)_(Mean)' 
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'NO'

    elif species == 'NO2':
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'NO2 mixing ratio (pptv)_(Mean)'
        unit_cut= 1e12 
        species_type = 'Conc.'
        actual_species_name = 'NO2'

    elif species == 'C2H6':
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'ethane mixing ratio (pptV)_(Mean)'
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'C2H6'

    elif species == 'C3H8':
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'propane mixing ratio (pptV)_(Mean)'
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'C3H8'

    elif species == 'DMS':
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'dms mixing ratio (pptV)_(Mean)'
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'DMS'

    elif species == 'TRA_6':  #Isoprene
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'Isoprene (pptv)_(Mean)'
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'Isoprene'

    elif species == 'ACET':
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'acetone mixing ratio (pptV)_(Mean)'
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'Acetone'

    elif species == 'GMAO_TEMP': # Temp from met fields
        units = 'K'
        first_label_pos = 3
        obs_data_name = 'Air Temperature (degC) Campbell_(Mean)'
        unit_cut= 1
        species_type = 'Temp.'
        actual_species_name = 'Surface Temperature'
        obs_switch = 1

    elif species == 'GMAO_PSFC': #Surface Pressure
        units = 'hPa'
        first_label_pos = 3
        obs_data_name = 'Atmospheric Pressure (hPa) Campbell_(Mean)'
        unit_cut= 1
        species_type = 'Pres.'
        actual_species_name = 'Surface Pressure'
        

    elif species == 'GMAO_WIND': #Wind Speed extirpolated from UWND and VWND 
        def read_diff_species():    
            k=names.index('GMAO_UWND')
            i=names.index('GMAO_VWND')
            model_cut=np.sqrt((model[:,k]**2)+(model[:,i]**2))
            return model_cut 
        units = r'$ms^{-1}$'
        first_label_pos = 3
        obs_data_name = 'Wind Speed (m/s) Campbell_(Mean)'
        unit_cut= 1
        species_type = 'Wind Speed'
        model_cut_switch = 1
        actual_species_name = 'Surface Windspeed'

    elif species == 'GMAO_RADSW': #Sensible heat flux form surface       
        units = r'$Wm^{-2}$'
        first_label_pos = 3
        obs_data_name = 'Solar Radiation (Wm-2) Campbell_(Mean)'
        unit_cut= 1
        species_type = 'Solar Radiation'
        actual_species_name = 'Surface Solar Radiation'

    elif species == 'GMAO_ABSH': #Absolute Humidity       
        units = 'molec/cm-3'
        first_label_pos = 3
        obs_data_name = ''
        unit_cut= 1
        species_type = 'Absolute Humidity'
        actual_species_name = 'Absolute Humidity'

    elif species == 'GMAO_RHUM': #Relative Humidity       
        units = '%'
        first_label_pos = 3
        obs_data_name = 'Relative Humidity (%) Campbell_(Mean)'
        unit_cut= 1
        species_type = 'Relative Humidity'
        actual_species_name = 'Relative Humidity'
    
   
#reads in the model data

    def readfile(filename, location):
        read = np.load(filename)
        names = read[0,2:]
        locs = np.where(read[:,1] == location)
        big = read[locs]
        valid_data = big[:,2:]
        names = names.tolist()
        valid_data = np.float64(valid_data)

        return valid_data, names

    try: 
        names
    except NameError:
# Readin the model output

        model , names = readfile("GEOS_v90103_4x5_CV_logs.npy","001") #001 represents CVO

# Processes the date 
        year=(model[:,0]//10000)
        month=((model[:,0]-year*10000)//100)
        day=(model[:,0]-year*10000-month*100)
    
        hour=model[:,1]//100
        min=(model[:,1]-hour*100)
    
        doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(2006,1,1,0,0,0) \
              for i in range(len(year))]
    
        since2006=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
        since2006=np.array(since2006)
       
        k=names.index(species)
        model_species = model[:,k]

        model_test = since2006 >= 730
        cut_time =   since2006[model_test]
        cut_model =  model_species[model_test]
       

        print cut_time 
        model_test_2 = cut_time < 1095 
        cut_time = cut_time[model_test_2]
        cut_model = cut_model[model_test_2]
        print len(cut_model)
        
        spacing =  1./24
        new_model_time = np.arange(0,365, spacing)
    
#now read in the observations

    myfile=nappy.openNAFile('York_merge_Cape_verde_1hr_R1.na')
    myfile.readData()

#ppy.openNAFile('York_merge_Cape_verde_1hr_R1.na')
    k_var1=myfile["VNAME"].index(obs_data_name)

# OK need to conver values from a list to a numpy array

    time=np.array(myfile['X'])
    
    if obs_switch == 0:
        var1=np.array(myfile['V'][k_var1])
    elif obs_switch == 1:
        var1=np.array(myfile['V'][k_var1])+273.15    

    valids1=var1 > 0

    time2=time[valids1]
     
    var2=var1[valids1]

    #cut time for 2007
 
    valids2 = time2 >= 365
    
    time2 = time2[valids2]
    var2 = var2[valids2]

    valids3 = time2 < 730
    time2 = time2[valids3]
    var2 = var2[valids3]
  
    print len(var2) 

    new_obs_time = [i-365 for i in time2]
           

#Pre normalise obs data for lomb analysis
    standard_deviation_obs_p = np.std(var2)
    mean_obs_p = np.mean(var2)
    normal_var2 = var2-mean_obs_p
    normal_var2 = normal_var2/standard_deviation_obs_p
 
#Calculate variance of pre-processed obs data- should be 1 if normal
    #standard_dev_obs = np.std(normal_var_2, dtype=np.float64)
    #variance_obs = standard_dev_obs**2
    #print 'Variance - pre-processed obs data= ', variance_obs

#Define sampling intervals
    samp_spacing = 1./24.

#Convert model time array into numpy array

#Need to normalise model data also
    if model_cut_switch == 0:
        k=names.index(species)
        model_cut = cut_model*unit_cut  	
    if model_cut_switch == 1:
        model_cut = read_diff_species()
    standard_deviation_model_p = np.std(model_cut)
    mean_model_p = np.mean(model_cut)
    normal_model = model_cut-mean_model_p
    normal_model = normal_model/standard_deviation_model_p

#Calculate variance of pre-processed model data- should be 1 if normal
    #standard_dev_model = np.std(normal_model, dtype=np.float64)
    #variance_model = standard_dev_model**2
    #print 'Variance - pre-processed model data= ', variance_model

#Plot them all up. 
    fig=plt.figure(figsize=(20,12))
    fig.patch.set_facecolor('white')
#Plot up standard conc. v time plot
    #ax1= fig.add_subplot(2, 1, 1)
    #fig.subplots_adjust(hspace=0.3)
    plt.plot(new_obs_time,var2, color='black', label='Cape Verde Obs. 2007')
    plt.plot(new_model_time, model_cut, color='green', label='GEOS v9.01.03 4x5 2008 ')
    plt.grid(True)
    leg=plt.legend(loc=first_label_pos, prop={'size':21})
    leg.get_frame().set_alpha(0.4)
    plt.xlabel('Time (Number of Days)', fontsize = 21)
    plt.ylim(0,70)
    plt.ylabel('%s (%s)' % (species_type,units), fontsize = 21)
    plt.title('%s V Time' % (actual_species_name), fontsize = 21)

#Define sampling frequency
    samp_freq = 24


    plt.show()



plot(species)
