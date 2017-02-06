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

species_list = ['CO'] #['CO','NO','NO2','C2H6','C3H8','DMS','TRA_6','ACET','GMAO_TEMP','GMAO_PSFC','GMAO_WIND','GMAO_RADSW','GMAO_RHUM']

Colour_list = ['green','black']



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

#Need to remove annual trend from data. Remove annual average from every value in dataset.
    #split into years

    months = [int(str(i)[4:6])  for i in model[:,0]]   
    months = np.array(months)
   
    month_count = 0
    present_month = months[0]
    monthly_array = []
    

    print len(model_cut)
   
    for i in model_cut:
        new_month = months[month_count]
        
        if present_month == new_month:
            monthly_array.append(i)
            present_month = new_month
            month_count+=1
        else: 
            monthly_ave =  np.average(monthly_array)
            detrended_month = [i + monthly_ave for i in monthly_array]
            month_total = 0
            monthly_array = []

            try:
                detrended_model_cut = np.hstack((detrended_model_cut,detrended_month))
            except:
                detrended_model_cut = np.array(detrended_month)
 
            monthly_array.append(i)
            present_month = new_month
            month_count+=1

        if month_count == len(model_cut):
            monthly_ave =  np.average(monthly_array)
            detrended_month = [i + monthly_ave for i in monthly_array]
            detrended_model_cut = np.hstack((detrended_model_cut,detrended_month))           

    print len(detrended_model_cut) 

    #test_jan = months == 01
    #test_feb = months == 02
    #test_mar = months == 03
    #test_apr = months == 04
    #test_may = months == 05
    #test_jun = months == 06
    #test_jul = months == 07
    #test_aug = months == 8
    #test_sep = months == 9
    #test_oct = months == 10
    #test_nov = months == 11
    #test_dec = months == 12

    #print test_jan

    #vals_jan = model_cut[test_jan]
    #vals_feb = model_cut[test_feb]
    #vals_mar = model_cut[test_mar]
    #vals_apr = model_cut[test_apr] 
    #vals_may = model_cut[test_may]
    #vals_jun = model_cut[test_jun]
    #vals_jul = model_cut[test_jul]
    #vals_aug = model_cut[test_aug]
    #vals_sep = model_cut[test_sep]
    #vals_oct = model_cut[test_oct]
    #vals_nov = model_cut[test_nov]     
    #vals_dec = model_cut[test_dec]

    #ave_jan = np.average(vals_jan)
    #ave_feb = np.average(vals_feb)
    #ave_mar = np.average(vals_mar)
    #ave_apr = np.average(vals_apr)
    #ave_may = np.average(vals_may)
    #ave_jun = np.average(vals_jun)
    #ave_jul = np.average(vals_jul)
    #ave_aug = np.average(vals_aug)
    #ave_sep = np.average(vals_sep)
    #ave_oct = np.average(vals_oct)
    #ave_nov = np.average(vals_nov)
    #ave_dec = np.average(vals_dec)

    #vals_jan = [i-ave_jan for i in vals_jan]
    #vals_feb = [i-ave_feb for i in vals_feb]
    #vals_mar = [i-ave_mar for i in vals_mar]
   #vals_apr = [i-ave_apr for i in vals_apr]
    #vals_may = [i-ave_may for i in vals_may]
    #vals_jun = [i-ave_jun for i in vals_jun]
    #vals_jul = [i-ave_jul for i in vals_jul]
    #vals_aug = [i-ave_aug for i in vals_aug]
    #vals_sep = [i-ave_sep for i in vals_sep]
    #vals_oct = [i-ave_oct for i in vals_oct]
    #vals_nov = [i-ave_nov for i in vals_nov]
    #vals_dec = [i-ave_dec for i in vals_dec] 

    #detrended_model_cut = np.hstack((vals_jan,vals_feb,vals_mar,vals_apr,vals_may,vals_jun,vals_jul,vals_aug,vals_sep,vals_oct,vals_nov,vals_dec)) 

    standard_deviation_model_p = np.std(model_cut)
    mean_model_p = np.mean(model_cut)
    normal_model = model_cut-mean_model_p
    normal_model = normal_model/standard_deviation_model_p

    standard_deviation_model_p = np.std(detrended_model_cut)
    mean_model_p = np.mean(detrended_model_cut)
    normal_model_2 = detrended_model_cut-mean_model_p
    normal_model_2 = normal_model_2/standard_deviation_model_p

    print normal_model_2    


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

#Model lomb
    detrend_fx, detrend_fy, nout, jmax, prob2 = lomb.fasper(since2006,normal_model_2, ofac, samp_freq)
#Divide output by sampling frequency
    detrend_fy = detrend_fy/samp_freq

#Create logged bins
    first_bin = []
    second_bin = []
    third_bin = []
    fourth_bin = []  
    fifth_bin = []
    sixth_bin = []

    freq_count = 0

    for i in detrend_fx:
        if i < 0.001:
            first_bin.append(detrend_fy[freq_count])
        if 0.001 <= i < 0.01:
            second_bin.append(detrend_fy[freq_count])
        if 0.01  <= i < 0.1:
            third_bin.append(detrend_fy[freq_count])
        if 0.1   <= i < 1:
            fourth_bin.append(detrend_fy[freq_count])
        if 1     <= i < 10:
            fifth_bin.append(detrend_fy[freq_count])
        if 10    <= i < 100:
            sixth_bin.append(detrend_fy[freq_count])  
        freq_count+=1

    print len(first_bin)
    print len(second_bin)
    print len(third_bin)
    print len(fourth_bin)
    print len(fifth_bin)
    print len(sixth_bin)

    if len(first_bin) < 10:
        second_bin = first_bin+second_bin
        first_bin = 0


    def bin_corrector(bin_name,next_bin_name,modulus):
        if modulus == 0:
            1+1
        else:
           selection = bin_name[-modulus:]
           bin_name = bin_name[:-modulus]
           next_bin_name =  selection + next_bin_name
        return bin_name, next_bin_name

# remainder of bins divided by 10 goes into previous bin, first and sixth order bins not full so do not treat same

    second_bin_mod = np.mod(len(second_bin),10)
    second_bin, third_bin = bin_corrector(second_bin, third_bin, second_bin_mod)
    third_bin_mod = np.mod(len(third_bin),10)
    third_bin, fourth_bin = bin_corrector(third_bin, fourth_bin, third_bin_mod)
    fourth_bin_mod =np.mod(len(fourth_bin),10)
    fourth_bin, fifth_bin = bin_corrector(fourth_bin, fifth_bin, fourth_bin_mod)
    fifth_bin_mod =np.mod(len(fifth_bin),10)
    fifth_bin, sixth_bin = bin_corrector(fifth_bin, sixth_bin, fifth_bin_mod)

    freq_count = 0
    second_freq = detrend_fx[:len(second_bin)]
    freq_count+=len(second_bin)
    third_freq = detrend_fx[freq_count:len(third_bin)]         
    freq_count+=len(third_bin)
    fourth_freq = detrend_fx[freq_count:len(fourth_bin)] 
    freq_count+=len(fourth_bin)
    fifth_freq = detrend_fx[freq_count:len(fifth_bin)] 
    freq_count+=len(fifth_bin)
    sixth_freq = detrend_fx[freq_count:]

    def chunks(l, n):
        return [l[i:i+n] for i in range(0, len(l), n)]


    if first_bin == 0:
        second_bin = chunks(second_bin,len(second_bin)/10)
        second_averages = np.average(second_bin, axis=1)
        second_freq = chunks(second_freq,len(second_freq)/10) 
        second_freq = np.average(second_freq, axis = 1)
        third_bin = chunks(third_bin,len(third_bin)/10)
        third_averages = np.average(third_bin, axis=1)
        third_freq = chunks(third_freq,len(third_freq)/10)
        third_freq = np.average(third_freq, axis = 1)
        fourth_bin = chunks(fourth_bin,len(fourth_bin)/10)
        fourth_averages = np.average(fourth_bin, axis=1)
        fourth_freq = chunks(fourth_freq,len(fourth_freq)/10)
        fourth_freq = np.average(fourth_freq, axis = 1) 
        fifth_bin_length = len(fifth_bin)
        fifth_bin = chunks(fifth_bin,len(fifth_bin)/10)
        fifth_averages = np.average(fifth_bin, axis=1)
        fifth_freq = chunks(fifth_freq,len(fifth_freq)/10)
        fifth_freq = np.average(fifth_freq, axis = 1)
  
        sixth_bin_n = (fifth_bin_length/10) * 10
        if sixth_bin_n > len(sixth_bin): 
            sixth_averages = np.average(sixth_bin)
            sixth_freq = np.average(sixth_freq)

        log_smoothed_powers =  np.hstack((second_averages,third_averages,fourth_averages,fifth_averages,sixth_averages))
        log_smoothed_freqs   = np.hstack((second_freq,third_freq,fourth_freq,fifth_freq,sixth_freq))  
     #else:
    

    #print len(first_bin)+len(second_bin)+len(third_bin)+len(fourth_bin)+len(fifth_bin)+len(sixth_bin)
 
    #smoothed_model_periods = model_periods[:len(smoothed_model)] 
    
    log_smoothed_periods = 1./log_smoothed_freqs

#plot up
    plt.loglog(1./fx, fy, color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='%s GEOS v9.01.03 4x5 ' %species)
    plt.loglog(log_smoothed_periods, log_smoothed_powers, color = 'green', marker = 'x', alpha = 0.75,markersize=2, label='GEOS 4x5 v90103 %s Smoothed ' %actual_species_name)
#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
    #nyquist_freq_lomb_model = frequencies[-1]
    #Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
    #print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 

#Obs. lomb
    #fa, fb, nout, jmax, prob2 = lomb.fasper(obs_time,normal_obs, ofac, samp_freq)
#Divide output by sampling frequency
    #fb = fb/samp_freq
    #reversed_fb = fb[::-1]

    #obs_periods = 1./fa
    #reversed_obs_periods = obs_periods[::-1]
    #cut_obs_periods, smoothed_obs = ema(reversed_fb,reversed_obs_periods,10,120)
    #smoothed_obs_periods = obs_periods[:len(smoothed_model)] 

#plot up
    #plt.loglog(1./fx, fy, color = colour_list[counter], marker = 'x', alpha = 0.75,markersize=2, label='%s GEOS v9.01.03 4x5 ' %species)
    #plt.loglog(cut_obs_periods, smoothed_obs, color = 'black', marker = 'x', alpha = 0.75,markersize=2, label='Cape Verde Obs. %s Smoothed ' %actual_species_name)

    plt.grid(True)
    leg=plt.legend(loc=4, prop={'size':21})
    leg.get_frame().set_alpha(0.4)
    plt.xlabel('Period (Days)')
    plt.ylabel(r'PSD $(ppb^{2}/days^{-1})$')
    plt.title('Lomb Scargle PSD.',fontsize=21)
    counter+=1



    plt.show()










