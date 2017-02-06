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

species_list = ['C2H6']

colour_list = ['green','black']



def readfile(filename, location):
    read = np.load(filename)
    names = read[0,2:]
    locs = np.where(read[:,1] == location)
    big = read[locs]
    valid_data = big[:,2:]
    names = names.tolist()
    valid_data = np.float64(valid_data)

    return valid_data, names

def obs_variable_finder(species):
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
    
 
    return units, obs_data_name, unit_cut, species_type, actual_species_name, obs_switch, model_cut_switch, ofac

def read_diff_species_wind():
    k=names.index('GMAO_UWND')
    i=names.index('GMAO_VWND')
    model_cut=np.sqrt((model[:,k]**2)+(model[:,i]**2))
    return model_cut

def ema(s, periods, n, period_limit):
    """
    returns an n period exponential moving average for
    the time series s

    s is a list ordered from oldest (index 0) to most
    recent (index -1)
    n is an integer

    returns a numeric array of the exponential
    moving average
    """
    s = np.array(s)
    ema = []
    cut_periods = []
    j = 1

    #get n sma first and calculate the next n period ema
    sma = sum(s[:n]) / n
    multiplier = 2 / float(1 + n)
    ema.append(sma)
    periods_mid = float(n/2)
    if np.mod(n,2) == 0:
        first_period = periods_mid-1
        valid_period = np.average(periods[first_period],periods[periods_mid])
        cut_periods.append(valid_period)
    else:
        valid_period = periods_mid-0.5
        cut_periods.append(periods[valid_period])     
 
    #EMA(current) = ( (Price(current) - EMA(prev) ) x Multiplier) + EMA(prev)
    ema.append(( (s[n] - sma) * multiplier ) + sma)
    cut_periods.append(periods[n])

    #now calculate the rest of the values
    periods_counter = n+1
    for i in s[n+1:]:
        if periods[periods_counter] < period_limit:
            #print 'below limit: ', i
            tmp = ( (i - ema[j]) * multiplier ) + ema[j]
            j = j + 1
            ema.append(tmp)
            cut_periods.append(periods[periods_counter])
            periods_counter +=1
            
    for i in periods[periods_counter:]:
        cut_periods.append(i)
    for i in s[periods_counter:]:
        ema.append(i)
            
    return cut_periods,ema

counter =0
#set plotting area. 
#fig=plt.figure(figsize=(20,12))
#fig.patch.set_facecolor('white')

try:
    names
except NameError:
# Readin the model output
    model , names = readfile("GEOS_logs.npy","001") #001 represents CVO
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

#Define sampling intervals
    samp_spacing = 1./24.

#Convert model time array into numpy array
    since2006=np.array(since2006)

for species in species_list:
    units, obs_data_name, unit_cut, species_type, actual_species_name, obs_switch, model_cut_switch, ofac = obs_variable_finder(species)

    fig=plt.figure(figsize=(20,12))
    fig.patch.set_facecolor('white')

#Need to normalise model data also
    if model_cut_switch == 0:
        k=names.index(species)
        model_cut = model[:,k]*unit_cut
    if model_cut_switch == 1:
        model_cut = read_diff_species_wind()

#Add daily emission trend onto ethane.

#Add dif. trend every min
#Add trend each hour

    first_quarter_day = np.linspace(0,30,num=6,endpoint=False)
    second_quarter_day = np.linspace(30,0,num=6,endpoint=False)
    third_quarter_day = np.linspace(0,-30,num=6,endpoint=False)
    fourth_quarter_day = np.linspace(-30,0,num=6,endpoint=False)
   
    daily_signal = np.concatenate((first_quarter_day,second_quarter_day,third_quarter_day,fourth_quarter_day)) 

    print daily_signal

    n=0
    n_end=24
   
    step=24 
    new_model_cut=[]
    hourly_count = 0

    model_data_length = len(model_cut)
    print model_data_length

    while n_end < model_data_length+1:
        
            sliced = model_cut[n:n_end]
            for a in sliced:
                new_model_cut.append(a+daily_signal[hourly_count])  
                hourly_count+=1
            n+=step
            n_end+=step
            hourly_count=0

    #print new_model_cut

    model_cut = new_model_cut
    print len(model_cut)

    standard_deviation_model_p = np.std(model_cut)
    mean_model_p = np.mean(model_cut)
    normal_model = model_cut-mean_model_p
    normal_model = normal_model/standard_deviation_model_p

#Calculate variance of pre-processed model data- should be 1 if normal
    #standard_dev_model = np.std(normal_model, dtype=np.float64)
    #variance_model = standard_dev_model**2
    #print 'Variance - pre-processed model data= ', variance_model

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

    obs_time=time[valids1]

    obs_var=var1[valids1]

#Pre normalise obs data for lomb analysis
    standard_deviation_obs_p = np.std(obs_var)
    mean_obs_p = np.mean(obs_var)
    normal_obs = obs_var-mean_obs_p
    normal_obs = normal_obs/standard_deviation_obs_p

#Calculate variance of pre-processed obs data- should be 1 if normal
    #standard_dev_obs = np.std(normal_var_2, dtype=np.float64)
    #variance_obs = standard_dev_obs**2
    #print 'Variance - pre-processed obs data= ', variance_obs


#Define sampling frequency
    samp_freq = 24

#Model lomb
    fx, fy, nout, jmax, prob2 = lomb.fasper(since2006,normal_model, ofac, samp_freq)
#Divide output by sampling frequency
    fy = fy/samp_freq
    reversed_fy = fy[::-1]  

 
    model_periods = 1./fx
    reversed_model_periods = model_periods[::-1]
    cut_model_periods, smoothed_model = ema(reversed_fy,reversed_model_periods,10,120) 
    #smoothed_model_periods = model_periods[:len(smoothed_model)] 
    
#plot up
    #plt.loglog(1./fx, fy, color = colour_list[counter], marker = 'x', alpha = 0.75,markersize=2, label='%s GEOS v9.01.03 4x5 ' %species)
    plt.loglog(cut_model_periods, smoothed_model, color = 'green', marker = 'x', alpha = 0.75,markersize=2, label='GEOS 4x5 v90103 %s Smoothed ' %actual_species_name)
#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
    #nyquist_freq_lomb_model = frequencies[-1]
    #Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
    #print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 

#Obs. lomb
    fa, fb, nout, jmax, prob2 = lomb.fasper(obs_time,normal_obs, ofac, samp_freq)
#Divide output by sampling frequency
    fb = fb/samp_freq
    reversed_fb = fb[::-1]

    obs_periods = 1./fa
    reversed_obs_periods = obs_periods[::-1]
    cut_obs_periods, smoothed_obs = ema(reversed_fb,reversed_obs_periods,10,120)
    #smoothed_obs_periods = obs_periods[:len(smoothed_model)] 

#plot up
    #plt.loglog(1./fx, fy, color = colour_list[counter], marker = 'x', alpha = 0.75,markersize=2, label='%s GEOS v9.01.03 4x5 ' %species)
    plt.loglog(cut_obs_periods, smoothed_obs, color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='Cape Verde Obs. %s Smoothed ' %actual_species_name)

    plt.grid(True)
    leg=plt.legend(loc=4, prop={'size':21})
    leg.get_frame().set_alpha(0.4)
    plt.xlabel('Period (Days)')
    plt.ylabel(r'PSD $(ppb^{2}/days^{-1})$')
    plt.title('Lomb Scargle PSD.',fontsize=21)
    counter+=1



    plt.show()










