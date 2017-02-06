import numpy as np
import csv
import glob
import datetime as datetime
import matplotlib.pyplot as plt
import logging as log
import scipy.stats as stats
import modules
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq, rfft
import collections

def readfile(filename, location):
    read = np.load(filename)
    names = read[0,2:]
    locs = np.where(read[:,1] == location)
    big = read[locs]
    valid_data = big[:,2:]
    names = names.tolist()
    valid_data = np.float64(valid_data)

    return valid_data, names

def readfile_GAW(filename, location):

	print location
	read = np.load(filename)
	names = read[0,1:]
	print names	
	locs = np.where(read[:,0] == location)
	big = read[locs]
	print big
	valid_data = big[:,1:]
	names = names.tolist()
	valid_data = np.float64(valid_data)

	return valid_data, names

def CV_data_reader(filename,obs_data_name, obs_switch):
	myfile=nappy.openNAFile(filename)
	myfile.readData()

	k_var1=myfile["VNAME"].index(obs_data_name)

# OK need to conver values from a list to a numpy array

	time=np.array(myfile['X'])

	if obs_switch == 0:
		var1=np.array(myfile['V'][k_var1])
	elif obs_switch == 1:
		var1=np.array(myfile['V'][k_var1])+273.15

	return time, var1
  
def NOAA_data_reader(filename):
	for files in filename:
		print files
		reader=csv.reader(open(files,'rb'), delimiter=' ', skipinitialspace=True)
		for row in reader:
			new = row[:]
			try:
				data.append(new)
			except:
				data=[new]

	data=np.array(data)
	year = data[:,1]
	month= data[:,2]
	day  = data[:,3]
	almost_date = [a+b for a,b in zip(year, month)]
	date = [a+b for a,b in zip(almost_date, day)]
	date = np.array(date)
	date = date.astype(int)
	hour = data[:,4]
	time = [i+'00' for i in hour]
	time = np.array(time)
	time = time.astype(int)
	vals = data[:,5]
	vals = np.float64(vals)
	return date, time, vals


def NOAA_data_reader_mace_head(filename):
    print filename
    reader=csv.reader(open(filename,'rb'), delimiter=',', skipinitialspace=True)
    for row in reader:
    	new = row[:]
        try:
            data.append(new)
        except:
            data=[new]

    data=np.array(data)
    date_p = data[:,0]
    year = [a[6:10] for a in date_p]
    month = [a[3:5] for a in date_p]
    day = [a[0:2] for a in date_p]
    almost_date = [a+b for a,b in zip(year, month)]
    date = [a+b for a,b in zip(almost_date, day)]
    date = np.array(date)
    date = date.astype(int)
    hour_p = data[:,1]
    hour = [i[0:2] for i in hour_p]
    counter =0
    for i in hour:
        if i == '24':
       	    date[counter] = date[counter+1]
            hour[counter] = '00'
        counter+=1
    time = [i+'00' for i in hour]
    time = np.array(time)
    time = time.astype(int)
    vals_p = data[:,2]
    vals=[]
    for i in vals_p:
        if i == '-':
            i = '-999'
    	vals.append(i)
    vals = np.array(vals)
    vals = vals.astype(np.float)
    #vals = vals.astype(int)
    return date, time, vals


def obs_variable_finder(species):
    model_cut_switch = 0
    obs_switch = 0
    ofac = 1
    if species == 'O3':
        units = 'ppbV'
        first_label_pos = 3
        obs_data_name = 'Ozone_mixing_ratio_ppbV_Mean'
        unit_cut= 1e9
        species_type = 'Conc.'
        actual_species_name = 'O3'
        ofac = 2.0001

    elif species == 'CO':
        units = 'ppbV'
        first_label_pos = 1
        obs_data_name = 'CO_mixing_ratio_ppbV_Mean'
        unit_cut= 1e9
        species_type = 'Conc.'
        actual_species_name = 'CO'
        ofac = 2.0001

    elif species == 'NO':
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'NO_mixing_ratio_pptv_Mean'
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'NO'

    elif species == 'NO2':
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'NO2_mixing_ratio_pptv_Mean'
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'NO2'

    elif species == 'C2H6':
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'ethane_mixing_ratio_pptV_Mean'
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'C2H6'

    elif species == 'C3H8':
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'propane_mixing_ratio_pptV_Mean'
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'C3H8'

    elif species == 'DMS':
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'dms_mixing_ratio_pptV_Mean'
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'DMS'

    elif species == 'TRA_6':  #Isoprene
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'Isoprene_pptv_Mean'
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'Isoprene'

    elif species == 'ACET':
        units = 'pptV'
        first_label_pos = 1
        obs_data_name = 'acetone_mixing_ratio_pptV_Mean'
        unit_cut= 1e12
        species_type = 'Conc.'
        actual_species_name = 'Acetone'

    elif species == 'GMAO_TEMP': # Temp from met fields
        units = 'K'
        first_label_pos = 3
        obs_data_name = 'Air_Temperature_degC_Campbell_Mean'
        unit_cut= 1
        species_type = 'Temp.'
        actual_species_name = 'Surface_Temperature'
        obs_switch = 1

    elif species == 'GMAO_PSFC': #Surface Pressure
        units = 'hPa'
        first_label_pos = 3
        obs_data_name = 'Atmospheric_Pressure_hPa_Campbell_Mean'
        unit_cut= 1
        species_type = 'Pres.'
        actual_species_name = 'Surface_Pressure'

    elif species == 'GMAO_WIND': #Wind Speed extirpolated from UWND and VWND 
        units = r'$ms^{-1}$'
        first_label_pos = 3
        obs_data_name = 'Wind_Speed_m/s_Campbell_Mean'
        unit_cut= 1
        species_type = 'Wind Speed'
        model_cut_switch = 1
        actual_species_name = 'Surface_Windspeed'

    elif species == 'GMAO_RADSW': #Sensible heat flux form surface       
        units = r'$Wm^{-2}$'
        first_label_pos = 3
        obs_data_name = 'Solar_Radiation_Wm-2_Campbell_Mean'
        unit_cut= 1
        species_type = 'Solar Radiation'
        actual_species_name = 'Surface_Solar_Radiation'

    elif species == 'GMAO_ABSH': #Absolute Humidity       
        units = 'molec/cm-3'
        first_label_pos = 3
        obs_data_name = ''
        unit_cut= 1
        species_type = 'Absolute Humidity'
        actual_species_name = 'Absolute_Humidity'

    elif species == 'GMAO_RHUM': #Relative Humidity       
        units = '%'
        first_label_pos = 3
        obs_data_name = 'Relative_Humidity_%_Campbell_Mean'
        unit_cut= 1
        species_type = 'Relative Humidity'
        actual_species_name = 'Relative_Humidity' 
    
 
    return units, obs_data_name, unit_cut, species_type, actual_species_name, obs_switch, model_cut_switch, ofac

def read_diff_species_wind(names,model):
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

def location_check(location):
    if location == 'Arrival_Heights':
        obsfile = 'obs_files/arrival_heights_o3_hourly/o3*'
        loc_label = 'Arrival Heights'
        model_index = '010'

    elif location == 'Barrow':
        obsfile = 'obs_files/barrow_o3_hourly/o3*'
        loc_label = 'Barrow'
        model_index = '015'

    elif location == 'Cape_Verde':
        obsfile = 'obs_files/York_merge_Cape_verde_1hr_R1.na' 
        loc_label = 'Cape Verde' 
        model_index = '036'  
 
    elif location == 'Lauder':
        obsfile = 'obs_files/lauder_o3_hourly/o3*'
        loc_label = 'Lauder'
        model_index = '106'

    elif location == 'Mace_Head':
        obsfile = 'obs_files/O3_mace_head_ppbV.txt'
        loc_label = 'Mace Head'
        model_index = '112'

    elif location == 'Mauna_Loa':
        obsfile = 'obs_files/mauna_loa_o3_hourly/o3*'
        loc_label = 'Mauna Loa'
        model_index = '116'

    elif location == 'Niwot_Ridge':
        obsfile = 'obs_files/niwot_ridge_o3_hourly/o3*'
        loc_label = 'Niwot Ridge'
        model_index = '132'

    elif location == 'Ragged_Point':
        obsfile = 'obs_files/ragged_point_o3_hourly/o3*'
        loc_label = 'Ragged Point'
        model_index = '152'

    elif location == 'South_Pole':
        obsfile = 'obs_files/south_pole_o3_hourly/o3*'
        loc_label = 'South Pole'
        model_index = '173'
 
    elif location == 'Trinidad_Head':
        obsfile = 'obs_files/trinidad_head_o3_hourly/o3*'
        loc_label = 'Trinidad Head'
        model_index = '189'

    elif location == 'Tudor_Hill':
        obsfile = 'obs_files/tudor_hill_o3_hourly/o3*'
        loc_label = 'Tudor Hill'
        model_index = '191'

    elif location == 'Tutuila':
        obsfile = 'obs_files/tutuila_o3_hourly/o3*'
        loc_label = 'Tutuila'
        model_index = '192'

    return obsfile,loc_label,model_index


def long_obs_data(location):
	if location == 'Mauna_Loa':
		obsfile = 'long_term_obs/Mauna_Loa/hourly_o3/o3*'
		loc_label = 'Mauna Loa'
	return obsfile, loc_label	
	
def date_process(date,time):
	year=(date//10000)
	month=((date-year*10000)//100)
	day=(date-year*10000-month*100)

	hour=time//100
	min=(time-hour*100)

	doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(2006,1,1,0,0,0) \
              for i in range(len(year))]

	processed_dates=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
	return processed_dates

def bin_corrector(bin_name,next_bin_name,modulus):
	if modulus == 0:
		1+1
	else:
		selection = bin_name[-modulus:]
		bin_name = bin_name[:-modulus]
	 	next_bin_name =  selection + next_bin_name
	return bin_name, next_bin_name

def chunks(l, n):
	return [l[i:i+n] for i in range(0, len(l), n)]


def binned_smooth(freq,power):
	first_bin = []
	second_bin = []
	third_bin = []
	fourth_bin = []  
	fifth_bin = []
	sixth_bin = []

	freq_count = 0

	for i in freq:
		if i < 0.001:
			first_bin.append(power[freq_count])
		if 0.001 <= i < 0.01:
			second_bin.append(power[freq_count])
		if 0.01  <= i < 0.1:
			third_bin.append(power[freq_count])
		if 0.1   <= i < 1:
			fourth_bin.append(power[freq_count])
		if 1     <= i < 10:
			fifth_bin.append(power[freq_count])
		if 10    <= i < 100:
			sixth_bin.append(power[freq_count])  
		freq_count+=1

	print len(first_bin)
	print len(second_bin)
	print len(third_bin)
	print len(fourth_bin)
	print len(fifth_bin)
	print len(sixth_bin)

#if len of first bin is less than 10, merge it with second bin.
	if len(first_bin) == 0:
		first_bin = 1	
		if len(second_bin) < 10:
			third_bin = second_bin+third_bin
			second_bin = 0
	elif len(first_bin) < 10:
		second_bin = first_bin+second_bin
		first_bin = 0

# remainder of bins divided by 10 goes into previous bin, first and sixth order bins not full so do not treat same
	if first_bin == 0:
		second_bin_mod = np.mod(len(second_bin),10)
		second_bin, third_bin = bin_corrector(second_bin, third_bin, second_bin_mod)
		third_bin_mod = np.mod(len(third_bin),10)
		third_bin, fourth_bin = bin_corrector(third_bin, fourth_bin, third_bin_mod)
		fourth_bin_mod =np.mod(len(fourth_bin),10)
		fourth_bin, fifth_bin = bin_corrector(fourth_bin, fifth_bin, fourth_bin_mod)
		fifth_bin_mod =np.mod(len(fifth_bin),10)
		fifth_bin, sixth_bin = bin_corrector(fifth_bin, sixth_bin, fifth_bin_mod)

		freq_count = 0
		second_freq = freq[:len(second_bin)]
		freq_count+=len(second_bin)
		third_freq = freq[freq_count:len(third_bin)]         
		freq_count+=len(third_bin)
		fourth_freq = freq[freq_count:len(fourth_bin)] 
		freq_count+=len(fourth_bin)
		fifth_freq = freq[freq_count:len(fifth_bin)] 
		freq_count+=len(fifth_bin)
		sixth_freq = freq[freq_count:]

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
	
	elif second_bin == 0:
		third_bin_mod = np.mod(len(third_bin),10)
		third_bin, fourth_bin = bin_corrector(third_bin, fourth_bin, third_bin_mod)
		fourth_bin_mod =np.mod(len(fourth_bin),10)
		fourth_bin, fifth_bin = bin_corrector(fourth_bin, fifth_bin, fourth_bin_mod)
		fifth_bin_mod =np.mod(len(fifth_bin),10)
		fifth_bin, sixth_bin = bin_corrector(fifth_bin, sixth_bin, fifth_bin_mod)

		freq_count = 0
		third_freq = freq[freq_count:len(third_bin)]
		freq_count+=len(third_bin)
		fourth_freq = freq[freq_count:len(fourth_bin)]
		freq_count+=len(fourth_bin)
		fifth_freq = freq[freq_count:len(fifth_bin)]
		freq_count+=len(fifth_bin)
		sixth_freq = freq[freq_count:]

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

		log_smoothed_powers =  np.hstack((third_averages,fourth_averages,fifth_averages,sixth_averages))
		log_smoothed_freqs   = np.hstack((third_freq,fourth_freq,fifth_freq,sixth_freq))

	else:
		print 'why3'
		first_bin_mod = np.mod(len(first_bin),10)
		first_bin, second_bin = bin_corrector(first_bin, second_bin, first_bin_mod)			
		second_bin_mod = np.mod(len(second_bin),10)
		second_bin, third_bin = bin_corrector(second_bin, third_bin, second_bin_mod)
		third_bin_mod = np.mod(len(third_bin),10)
		third_bin, fourth_bin = bin_corrector(third_bin, fourth_bin, third_bin_mod)
		fourth_bin_mod =np.mod(len(fourth_bin),10)
		fourth_bin, fifth_bin = bin_corrector(fourth_bin, fifth_bin, fourth_bin_mod)
		fifth_bin_mod =np.mod(len(fifth_bin),10)
		fifth_bin, sixth_bin = bin_corrector(fifth_bin, sixth_bin, fifth_bin_mod)

		freq_count = 0
		first_freq = freq[:len(first_bin)]
		freq_count+=len(first_bin)
		second_freq = freq[:len(second_bin)]
		freq_count+=len(second_bin)
		third_freq = freq[freq_count:len(third_bin)]         
		freq_count+=len(third_bin)
		fourth_freq = freq[freq_count:len(fourth_bin)] 
		freq_count+=len(fourth_bin)
		fifth_freq = freq[freq_count:len(fifth_bin)] 
		freq_count+=len(fifth_bin)
		sixth_freq = freq[freq_count:]

		first_bin = chunks(first_bin,len(first_bin)/10)
		first_averages = np.average(first_bin, axis=1)
		first_freq = chunks(first_freq,len(first_freq)/10)
		first_freq = np.average(first_freq, axis = 1)
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

		log_smoothed_powers =  np.hstack((first_averages,second_averages,third_averages,fourth_averages,fifth_averages,sixth_averages))
		log_smoothed_freqs   = np.hstack((first_freq,second_freq,third_freq,fourth_freq,fifth_freq,sixth_freq))  

	log_smoothed_periods = 1./log_smoothed_freqs

	return log_smoothed_periods, log_smoothed_powers


def ema_keep_peaks(s, periods, n, period_limit):
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
	multiplier = 2 / float (1 + n)
	ema.append(sma)
	periods_mid = float(n/2)
	
	if np.mod(n,2) == 0:
		periods_mid = int(periods_mid)
		period_mid_val =  periods[periods_mid]
		first_period = periods_mid-1
		first_period_val =  periods[first_period]
		valid_period = np.average((first_period_val,period_mid_val))
		cut_periods.append(valid_period)
	else:
		valid_period = periods_mid-0.5
		cut_periods.append(periods[valid_period])

	#EMA(current) = ( (Price(current) - EMA(prev) ) x Multiplier) + EMA(prev)
	ema.append(( (s[n] - sma) * multiplier ) + sma)
	cut_periods.append(periods[n])

	adjust_f = (n-len(ema))+1


	daily_period = 1
	half_annual_period = (365.25/2)
	annual_period = 365.25

	closest_daily_period_index = min(range(len(periods)), key=lambda i: abs(periods[i]-daily_period))
	closest_half_annual_period_index = min(range(len(periods)), key=lambda i: abs(periods[i]-half_annual_period))
	closest_annual_period_index = min(range(len(periods)), key=lambda i: abs(periods[i]-annual_period))

	daily_lower_index, daily_upper_index = key_period_check(s,closest_daily_period_index)
	half_annual_lower_index, half_annual_upper_index = key_period_check(s,closest_half_annual_period_index)
	annual_lower_index, annual_upper_index = key_period_check(s,closest_annual_period_index)

	
	daily_array = []
	half_annual_array = []
	annual_array = []
	daily_periods = []
	half_annual_periods = []
	annual_periods = []		

	daily_count = 0
	half_annual_count = 0
	annual_count = 0

	#now calculate the rest of the values
	periods_counter = n+1
	for i in s[n+1:]:
		if periods[periods_counter] < period_limit:
			current_period_counter = periods_counter
			if (periods_counter >= daily_lower_index) & (periods_counter <= daily_upper_index):
				tmp = i
				#j = j + 1
				#cut_periods.append(periods[periods_counter])
				#ema.append(tmp)
				daily_array.append(tmp)
				daily_periods.append(periods[periods_counter])
				if daily_count == 0:
					daily_index = periods_counter
				periods_counter +=1
				daily_count+=1
			elif (periods_counter >= half_annual_lower_index) & (periods_counter <= half_annual_upper_index):
				tmp = i
				#j = j + 1
				#cut_periods.append(periods[periods_counter])
				#ema.append(tmp)
				half_annual_array.append(tmp)
				half_annual_periods.append(periods[periods_counter])
				if half_annual_count == 0:
					half_annual_index = periods_counter
				periods_counter +=1
				half_annual_count+=1
			elif (periods_counter >= annual_lower_index) & (periods_counter <= annual_upper_index):
				tmp = i
				#j = j + 1
				#cut_periods.append(periods[periods_counter])
				#ema.append(tmp)
				annual_periods.append(periods[periods_counter])
				annual_array.append(tmp)
				if annual_count == 0:
					annual_index = periods_counter
				periods_counter +=1
				annual_count+=1
			else:
				tmp = ( (i - ema[j]) * multiplier ) + ema[j]
				j = j + 1
				ema.append(tmp)
				cut_periods.append(periods[periods_counter])
				periods_counter +=1

	for i in periods[periods_counter:]:
		cut_periods.append(i)
	for i in s[periods_counter:]:
		ema.append(i)

	#reinsert key periods back into arrays
	#daily
	cut_periods[daily_index-adjust_f:daily_index-adjust_f] =  daily_periods
	ema[daily_index-adjust_f:daily_index-adjust_f] = daily_array
	#half_annual
	cut_periods[half_annual_index-adjust_f:half_annual_index-adjust_f] = half_annual_periods
	ema[half_annual_index-adjust_f:half_annual_index-adjust_f] = half_annual_array	
	#annual
	cut_periods[annual_index-adjust_f:annual_index-adjust_f] = annual_periods
	ema[annual_index-adjust_f:annual_index-adjust_f] = annual_array	


	return cut_periods,ema

def key_period_check(input_list,closest_index):
 #function to test and leave significant peaks
	test_downwards_index = 0
	test_upwards_index = 0

	#vals = input_list[closest_index-test_downwards_index:closest_index+test_upwards_index]
	vals = input_list[closest_index]
	mean = np.average(vals)
	percent_range = 90
	mean_percent = (mean/100) * percent_range
	mean_percent = np.abs(mean_percent)
	bottom_index = closest_index-test_downwards_index
	top_index = closest_index+test_upwards_index
    #start to check if values outside first 3 indices are 10% outside mean - if so the check stops, if not it continues
	manip_vals = vals
	prev_val = vals

	test = 0
    #test downwards
	while test==0:
		test_downwards_index +=1	
		testdown_val = input_list[closest_index-test_downwards_index]
		if (testdown_val > (mean-mean_percent)) & (testdown_val < prev_val) :      
			manip_vals = np.insert(manip_vals,0,testdown_val)
			#mean = np.average(manip_vals)
			#mean_percent = (mean/100) * percent_range

			bottom_index = closest_index-test_downwards_index
		else:
			test+=1
		
		prev_val = input_list[closest_index-test_downwards_index]

	test = 0	
	manip_vals = vals
	mean = np.average(manip_vals)
	prev_val = vals

    #test upwards 
	while test==0:
		test_upwards_index  +=1
		try:
			testup_val = input_list[closest_index+test_upwards_index]
             
			if (testup_val > (mean-mean_percent)) & (testup_val < prev_val):
				manip_vals = np.insert(manip_vals,-1,testup_val)
			#mean = np.average(manip_vals)
			#mean_percent = (mean/100) * percent_range

				top_index = closest_index+(test_upwards_index+1)	
			else:
				test+=1

			prev_val = input_list[closest_index+test_upwards_index]
		except:
			test+=1


	return bottom_index, top_index


def significant_peaks(x,y,running_window_len=20,ahead=30,sig = 2.5):

	counter = 0
	window_start = 1
	
	sig_indices = []
	sig_periods = []
	sig_mag = []	

	
	for current_mag in y:
		#print counter
		if counter == 0:
			running_window = current_mag
			counter+=1		
			#print '1'	
		elif (counter > 0) & (counter < running_window_len):
			#test if current_mag is certain factor greater than running_window
			if current_mag > running_window*sig:
				#test ahead to see if actual peak or just a jitter
				ahead_array = y[counter+1:counter+(1+ahead)] 
				ahead_array_ave = np.average(ahead_array)				
				#if current_mag > ahead_array_ave*sig:
				sig_periods.append(x[counter])
				sig_mag.append(current_mag)
				sig_indices.append(counter)
			# if peak is identified, don't include it in running_window 
				running_window =  running_window
				counter+=1
				#else:   
				#	running_window = np.average([running_window,current_mag])
				#	counter+=1 
			else:	
				running_window = np.average([running_window,current_mag])
				counter+=1		

			#print '2'

		elif (counter >= running_window_len):
			select_window=[]
			remove_indices = []
			test_array = []

			#if x[counter] > 300:
			#	print 'period = ', x[counter]
			#	print 'current mag = ',current_mag
			#	print 'limit= ', running_window*sig
			#test if current_mag is certain factor greater than running_window
			if current_mag > running_window*sig:
				#test ahead to see if actual peak or just a jitter
				ahead_array = y[counter+1:counter+(1+ahead)]          
				ahead_array_ave = np.average(ahead_array)
				#if current_mag > ahead_array_ave*sig:
				sig_periods.append(x[counter])
				sig_mag.append(current_mag)
				sig_indices.append(counter)
				
				#remove magnitudes identified as peak from window	
				running_window = running_window
				counter+=1			
				window_start+=1

				#else:
				#	test_array = np.arange(window_start,window_start+running_window_len,1)
				#	test_array = y[test_array]					
				#	running_window = np.average(test_array)
				#	counter+=1
				#	window_start+=1

			else:
				test_array = np.arange(window_start,window_start+running_window_len,1)
				test_array = y[test_array]
				running_window = np.average(test_array)
				counter+=1
				window_start+=1	
	
	#remove duplicate periods

	previous_period = 0
	del_indices = []
	for i in range(len(sig_periods)):
		if sig_periods[i] < 0.95:
			del_indices.append(i)
		else:
			period_percent = sig_periods[i] / 10
			if (previous_period > (sig_periods[i]-period_percent)) & (previous_period < (sig_periods[i]+period_percent)):
				if previous_mag >= sig_mag[i]:
					del_indices.append(i)
				elif sig_mag[i] > previous_mag:   
					del_indices.append(i-1)
		
		previous_period = sig_periods[i]
		previous_mag = sig_mag[i] 		
	
	sig_periods = np.delete(sig_periods,del_indices)			
	sig_mag = np.delete(sig_mag,del_indices)
	sig_indices = np.delete(sig_indices,del_indices)

	ns = 0	
	nt = 1
	tester = True
	temp_array = [ns]
	del_indices2 = []
	
	for a in range(len(sig_periods)-1):
		percentile = sig_periods[ns]/10
		nt = 1+ns
		while tester == True:
			if (sig_periods[nt] - percentile) < sig_periods[ns]:
				temp_array.append(nt)
				nt+=1
				ns+=1
				percentile = sig_periods[ns]/10
			else:
				tester = False
		if len(temp_array) > 1:
			tester = True
			ns= a+1
			#print 'periods = ', sig_periods[temp_array]
			#print 'mags = ', sig_mag[temp_array]
			temp_mags = sig_mag[temp_array]
			temp_mag_max_index = max(enumerate(temp_mags),key=lambda x: x[1])[0]
			temp_mag_max_index= temp_array[temp_mag_max_index]
			#print 'max index = ', temp_mag_max_index
			#print 'temp array = ', temp_array
			for x in temp_array:
				if x != temp_mag_max_index:
					del_indices2.append(x)	
					#print 'x = ', x				
			#print 'del indices = ',  del_indices2
			temp_array = [ns]
		else:
			tester = True			
			ns= a+1
			temp_array = [ns]

	del_indices2 = np.unique(del_indices2)
	sig_periods = np.delete(sig_periods,del_indices2)
	sig_mag = np.delete(sig_mag,del_indices2)
	sig_indices = np.delete(sig_indices, del_indices2)

	print 'Significant Peak Periods = ', sig_periods
	return sig_periods, sig_mag, sig_indices

def align_periods(obs_time,model_time,obs_vals,model_cut):
#cut obs and model data to be for the same times
	#check if obs. start time is same as model start time
	if obs_time[0] != model_time[0]:
		#determine which time is higher
		if obs_time[0] > model_time[0]:
			start_ref = 'obs'
		if model_time[0] > obs_time[0]:
			start_ref = 'model'

		if start_ref == 'obs':
			start_point = obs_time[0]
			new_start_point = np.argmin(np.abs(model_time - start_point))
			model_time = model_time[new_start_point:]
			model_cut = model_cut[new_start_point:]
		else:
			start_point = model_time[0]
			new_start_point = np.argmin(np.abs(obs_time - start_point))
			obs_time = obs_time[new_start_point:]
			obs_vals = obs_vals[new_start_point:]


	#check if obs. end time is same as model end time
	if obs_time[-1] != model_time[-1]:
	#determine which time is lower
		if obs_time[-1] > model_time[-1]:
			end_ref = 'model'
		if model_time[-1] > obs_time[-1]:
			end_ref = 'obs'

		if end_ref == 'obs':
			end_point = obs_time[-1]
			new_end_point = np.argmin(np.abs(model_time - end_point))
			model_time = model_time[:new_end_point+1]
			model_cut = model_cut[:new_end_point+1]
		else:
			end_point = model_time[-1]
			new_end_point = np.argmin(np.abs(obs_time - end_point))
			obs_time = obs_time[:new_end_point+1]
			obs_vals = obs_vals[:new_end_point+1]

	return obs_time,model_time,obs_vals,model_cut

def read_CV(filename,obs_data_name, obs_switch):
	data = np.genfromtxt(filename,skip_header = 428, names=True)
	names = data.dtype.names
	time = data['Day_of_yearstart']
	if obs_switch == 0:
		var1=data[obs_data_name]
	elif obs_switch == 1:
		var1=data[obs_data_name]+273.15

	return time, var1

def automatic_smooth(s, periods, n, sig_indices):
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
    multiplier = 2 / float (1 + n)
    ema.append(sma)
    periods_mid = float(n/2)

    if np.mod(n,2) == 0:
        periods_mid = int(periods_mid)
        period_mid_val =  periods[periods_mid]
        first_period = periods_mid-1
        first_period_val =  periods[first_period]
        valid_period = np.average((first_period_val,period_mid_val))
        cut_periods.append(valid_period)
    else:
        valid_period = periods_mid-0.5
        cut_periods.append(periods[valid_period])

    #EMA(current) = ( (Price(current) - EMA(prev) ) x Multiplier) + EMA(prev)
    ema.append(( (s[n] - sma) * multiplier ) + sma)
    cut_periods.append(periods[n])

    #start_index_list = np.zeros(len(sig_indices))
    #end_index_list = np.zeros(len(sig_indices))

    #point = 0
    #if len(sig_indices) != 0:
        #for index in sig_indices:
            #start_index_list[point], end_index_list[point] = key_period_check(s,index)
            #point+=1

	#test if any indices are repeated.
    #diff_array = [b-a for a,b in zip(start_index_list,end_index_list)]

    #append_indices = []   
 
    #for a in range(len(end_index_list)):
     #   appender = []
      #  diff = diff_array[a]
       # if diff == 0:
        #    diff = 1
        #for b in np.arange(diff):
         #   end_i = end_index_list[a]
          #  toappend = end_i - b
           # appender.append(toappend)
        #append_indices.append(appender)
    
    
    #first = 0
    #second = 0
    #remove_indices = []

    #new_append_indices = []


    #for i in append_indices:
     #   if len(i) > 1:
      #      new_append_indices.append(i[1:])
		
    #append_indices = new_append_indices

    #for i in range(len(append_indices)-1):
       	#first_array = append_indices[first]
        #print append_indices[first]
       # second = i+1
       # for x in range(len(append_indices)-second):
          #  second_array = append_indices[second]
         #   matches = set(first_array) & set(second_array)
          #  n_matches = len(matches)
          #  if n_matches > 0:
             #   if np.max(first_array) > np.max(second_array):
               #     remove_indices.append(second)
              #  elif np.max(first_array) == np.max(second_array):
              #      a = np.sum(first_array)
             #       b = np.sum(second_array)
            #        if a >= b:
           #             remove_indices.append(second)
          #          else:
         #               remove_indices.append(first)
        #        else:
       #             remove_indices.append(first)
      #      second+=1 
     #   first+=1


    #start_index_list = np.delete(start_index_list,remove_indices)
    #end_index_list = np.delete(end_index_list,remove_indices)



    #start_index_list = start_index_list.astype(int) 
    #end_index_list = end_index_list.astype(int)
    #print 'start periods', periods[start_index_list][0:-1]	
    #print 'end periods', periods[end_index_list][0:-1]



# take square root of addition of squares in peak periods
#    for n in range(len(start_index_list)):
#        if end_index_list[n] - start_index_list[n] > 1:
 #           temp_arr = s[start_index_list[n]:end_index_list[n]]
  #          temp_index_arr = np.arange(start_index_list[n],end_index_list[n])
   #         arr_max = max(temp_arr)
    #        temp_index = np.where(temp_arr==arr_max)
     #       max_index = temp_index_arr[temp_index]
            
      #      negative_test = temp_arr < 0 
       #     square_vals = np.square(temp_arr)
        #    square_vals[negative_test] = -np.abs(square_vals[negative_test])        

         #   square_sum = np.sum(square_vals)
          #  if square_sum < 0:
           #     temp_sum = np.abs(square_sum)
            #    temp_sqrt = np.sqrt(temp_sum)
             #   sqrt_sum = -np.abs(temp_sqrt)
           # else:
            #    sqrt_sum = np.sqrt(square_sum)

            #start_index_list[n] = max_index
           # end_index_list[n] = max_index 

            #max_index = max_index[0]

            #s[max_index] = sqrt_sum  


    #full_index_list= []

    #for a,b in zip(start_index_list,end_index_list):
     #   full_index_list.append(a)
      #  full_index_list.append(b)
    #print 'full index list', full_index_list 

    #if len(sig_indices) != 0:
     #   bot = 0
      #  top = 1
       # current_bottom_index = full_index_list[bot]
        #current_top_index = full_index_list[top]

    #elif len(sig_indices) == 0:
     #   current_bottom_index = -99	
      #  current_top_index = -99        

    append_array = []
    append_periods = []

    key_index = []
    end_index = []

    s_index = []
    e_index = []    

    p_count = 0

    s_ref = 0
    e_ref = 1

#now calculate the rest of the values
    periods_counter = n+1
    for i in s[n+1:]:
        if periods[periods_counter] < period_limit:
            current_period_counter = periods_counter
			
            if (periods_counter >= current_bottom_index) & (periods_counter <= current_top_index):
                tmp = i
                append_array.append(tmp)
                append_periods.append(periods[periods_counter])
                if p_count == 0:
                    key_index.append(periods_counter)
                    s_index.append(s_ref)
                periods_counter +=1
                p_count+=1
                if periods_counter > current_top_index:
                    e_index.append(e_ref)
                    s_ref = e_ref
                    e_ref+=1
                    try:
                        bot+=2
                        top+=2
                        current_bottom_index = full_index_list[bot]
                        current_top_index = full_index_list[top]
                        end_index.append(periods_counter)
                        p_count=0
                    except:
                        1+1
                else:
                    e_ref+=1	
				

            else:
                tmp = ( (i - ema[j]) * multiplier ) + ema[j]
                j = j + 1
                ema.append(tmp)
                cut_periods.append(periods[periods_counter])
                periods_counter +=1
							
    for i in periods[periods_counter:]:
        cut_periods.append(i)
    for i in s[periods_counter:]:
        ema.append(i)
	
    #reinsert key periods back into arrays
    if len(sig_indices) != 0:
        for num in range(len(key_index)):
            cut_periods = np.append(cut_periods,append_periods[s_index[num]])
            ema = np.append(ema,append_array[s_index[num]])
        sorted_period_indices = sorted(range(len(cut_periods)), key=lambda k: cut_periods[k])
        cut_periods = cut_periods[sorted_period_indices]
        ema = ema[sorted_period_indices]

    return cut_periods,ema

def peak_rank(x,y):
	
	#amend_y = y

	#counter = 0
	#window_start = 1

	#for current_mag in y:
	#	if counter == 0:
	#		running_window = current_mag
	#		counter+=1
	#	elif (counter > 0) & (counter < running_window_len):
	#		running_window = np.average([running_window,current_mag])
	#		counter+=1
	#		amend_y[counter] = y[counter] - running_window

	#	elif (counter >= running_window_len):
            #remove_indices = []
	#		test_array = []

	#		test_array = np.arange(window_start,window_start+running_window_len,1)
	#		test_array = y[test_array]
	#		running_window = np.average(test_array)
	#		counter+=1
	#		window_start+=1
	#		if counter == len(y):
	#			1+1
	#		else:
	#			amend_y[counter] = y[counter] - running_window


	#print amend_y
	#print len(amend_y)
	#print len(y)
#sort arrays into descending order by magnitude
	sorted_y_indices = sorted(range(len(y)), key=lambda k: y[k],reverse=True)
	sorted_y = y[sorted_y_indices]
	sorted_x = x[sorted_y_indices]

	top_y = sorted_y[0:50]
	top_x = sorted_x[0:50]
	top_indices = sorted_y_indices[0:50]
	top_y = np.array(top_y)
	top_x = np.array(top_x)
	top_indices = np.array(top_indices)
	print 'top y ', top_y
	print 'top x ', top_x
	print 'top indices ', top_indices

	start = 1
	found = True
	test_n = 0
	full_count=0
	remove_indices = []

#remove duplicate periods
	for indi in top_indices:
		print 'indi = ',indi
		mini_array = []	
		while found == True:
			n_found = 0
			counter = start
			test_n+=1
			for indi_next in top_indices[start:]:
				#print 'indi', indi
				#print 'indi_next',indi_next
				#print 'test_n',test_n
				if (indi_next == indi-test_n) or (indi_next == indi+test_n):
					#print indi,indi_next
					mini_array.append(counter)
					#print counter
					n_found+=1
				counter+=1
			if n_found == 0:
				print 'should exit'
				found = False
		if len(mini_array) > 0:
			mini_array.append(full_count)
			mini_array_vals = top_y[mini_array]
			max_index = np.argmax(mini_array_vals)
			max_index = mini_array[max_index]
			#print top_indices[mini_array]
			print mini_array
			print max_index
			for ii in mini_array:
				if ii != max_index:
					remove_indices.append(ii)
	
		start+=1
		full_count+=1
		found = True
		test_n = 0			

	remove_indices = list(set(remove_indices))
	top_y = np.delete(top_y,remove_indices)
	top_x = np.delete(top_x,remove_indices)	
	top_indices = np.delete(top_indices,remove_indices)

	print 'top y ', top_y
	print 'top x ', top_x
	print 'top indices ', top_indices

#remove key peaks not within certain range of top peak
	top_peak_index = np.argmax(top_y)
	top_peak = top_y[top_peak_index]
	limit = (top_peak/100)*50 #peaks must be with 75% of top peak or else removed
	limit = top_peak - limit
	remove_peaks = []
	for counter in range(len(top_y)):
		if top_y[counter] < limit:
			remove_peaks.append(counter)

	top_y = np.delete(top_y,remove_peaks)
	top_x = np.delete(top_x,remove_peaks)
	top_indices = np.delete(top_indices,remove_peaks)

	#print remove_indices

	print 'top y ', top_y
	print 'top x ', top_x
	print 'top indices ', top_indices

#start+=1
# sort cut arrays into ascending order by periods
	re_sorted_x_indices = sorted(range(len(top_x)), key=lambda k: top_x[k])
	sig_x = top_x[re_sorted_x_indices]
	sig_y = top_y[re_sorted_x_indices]
	sig_indices = top_indices[re_sorted_x_indices]

	print 'Significant Periods = ', sig_x
	return sig_x, sig_y, sig_indices
	
def automatic_smooth_2(s, periods, n, sig_indices):
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
	multiplier = 2 / float (1 + n)
	ema.append(sma)
	periods_mid = float(n/2)

	if np.mod(n,2) == 0:
		periods_mid = int(periods_mid)
		period_mid_val =  periods[periods_mid]
		first_period = periods_mid-1
		first_period_val =  periods[first_period]
		valid_period = np.average((first_period_val,period_mid_val))
		cut_periods.append(valid_period)
	else:
		valid_period = periods_mid-0.5
		cut_periods.append(periods[valid_period])

    #EMA(current) = ( (Price(current) - EMA(prev) ) x Multiplier) + EMA(prev)
	ema.append(( (s[n] - sma) * multiplier ) + sma)
	cut_periods.append(periods[n])

	append_array = []
	append_periods = []
	key_index = []
	s_index = []

#now calculate the rest of the values
	sig_index = 0
	key_period_index = sig_indices[sig_index]

	periods_counter = n+1
	for i in s[n+1:]:
		current_period_counter = periods_counter

		if periods_counter == key_period_index:
			tmp = i
			append_array.append(tmp)
			append_periods.append(periods[periods_counter])
			key_index.append(periods_counter)
			periods_counter +=1
			try:
				sig_index+=1
				key_period_index = sig_indices[sig_index]
			except:
				1+1

		else:
			tmp = ( (i - ema[j]) * multiplier ) + ema[j]
			j = j + 1
			ema.append(tmp)
			cut_periods.append(periods[periods_counter])
			periods_counter +=1

	for i in periods[periods_counter:]:
		cut_periods.append(i)
	for i in s[periods_counter:]:
		ema.append(i)	

    #reinsert key periods back into arrays
	if len(sig_indices) != 0:
		for num in range(len(key_index)):
			cut_periods = np.append(cut_periods,append_periods[num])
			ema = np.append(ema,append_array[num])
		sorted_period_indices = sorted(range(len(cut_periods)), key=lambda k: cut_periods[k])
		cut_periods = cut_periods[sorted_period_indices]
		ema = ema[sorted_period_indices]
	

	return cut_periods,ema
	
def grid_obs_centre_convergance(lat_e,lon_e,obs_lats,obs_lons): 
#check which gridboxes each obs_lat & obs_lon convergance point lie in & then return central convergance point of lat and lon in that gridbox
	obs_lats_centre = []
	obs_lons_centre = []
	model_indices = []

	for a,b in zip(obs_lats, obs_lons):
		lat_i = np.searchsorted(lat_e,a,side='left')
		lon_i = np.searchsorted(lon_e,b,side='left')

		lower_lat = lat_e[lat_i-1]
		upper_lat = lat_e[lat_i]
		lower_lon = lon_e[lon_i-1]
		upper_lon = lon_e[lon_i]

		lat_diff = upper_lat - lower_lat
		lon_diff = upper_lon - lower_lon
		centre_lat = lower_lat + (lat_diff/2.)
		centre_lon = lower_lon + (lon_diff/2.)

		obs_lats_centre=np.append(obs_lats_centre,centre_lat)
		obs_lons_centre=np.append(obs_lons_centre,centre_lon)
		model_indices = np.append(model_indices,[int(lat_i)-1,int(lon_i)-1])
	model_indices = np.reshape(model_indices, (-1,2))
	model_indices = model_indices.astype(int)
	
	return obs_lats_centre, obs_lons_centre, model_indices 
	
def correct_grid_daily_phase(model,lon_c,lat_c,lon_e,lat_e):

    start = 0
    end = len(lon_c)

    adjusted_values = []
    correction = (np.pi*2)/len(lon_c)

    first_array = np.arange(len(lon_c)/2,len(lon_c),1)
    second_array = np.arange(0,len(lon_c)/2,1)

    time_correction = np.append(first_array,second_array)

    normal_time = range(len(lon_c))

	#correct model phase

    for num in range(len(time_correction)):
        longitude_band = model[normal_time[num]::len(lon_c)]

        for i in longitude_band:
            val = i-(time_correction[num]*correction)
            if val < (-np.pi):
                val = np.abs(val)
                remain = val - np.pi
                val = np.pi - remain
                adjusted_values = np.append(adjusted_values,val)
            elif val > (np.pi):
                val = np.abs(val)
                remain = val - np.pi
                val = -np.pi + remain
                adjusted_values = np.append(adjusted_values,val)
            else:
                adjusted_values = np.append(adjusted_values,val)

	#reshape model vals into latitude bands
    adjusted_values_2 = []

    for start in range(len(lat_c)):
        latitude_band = adjusted_values[start::len(lat_c)]	
        adjusted_values_2 = np.append(adjusted_values_2,latitude_band)

    model = adjusted_values_2
    return model

def correct_select_daily_phase(obs,lon_c,lat_c,lon_e,lat_e, obs_lons):
    correction = (np.pi*2)/len(lon_c)

#correct obs daily phase output
#for each value, produce a correction number depending on which longitude box the value is in and then correct  
    for num in range(len(obs)):
        lon_i = np.searchsorted(lon_e,obs_lons[num],side='left')
        lon_i-=1
        if lon_i <len(lon_c)/2:
            lon_i = lon_i+len(lon_c)/2
        elif lon_i >=len(lon_c)/2:
            lon_i = lon_i-len(lon_c)/2
        obs[num] = obs[num]-(lon_i*correction)
        if obs[num] < (-np.pi):
            val = np.abs(obs[num])
            remain = val - np.pi
            val = np.pi - remain
            obs[num] = val
        elif obs[num] > (np.pi):
            val = np.abs(obs[num])
            remain = val - np.pi
            val = -np.pi + remain
            obs[num] = val
    return obs
	
def convert_phase_units(obs,model,con_num):
	for num in range(len(model)):
		if model[num] < 0:
			model[num] = np.pi- np.abs(model[num])
		elif model[num] > 0:
			model[num] = np.pi + model[num]
				
	for num in range(len(obs)):
		if obs[num] < 0:
			obs[num] = np.pi- np.abs(obs[num])
		elif obs[num] > 0:
			obs[num] = np.pi + obs[num]

	convert_factor = con_num/(2*np.pi)
	model = model*convert_factor
	obs = obs*convert_factor
	
	return obs, model

def grid_reshape(lon_c,lat_c,model):
    start = 0
    end = len(lon_c)

    for i in range(len(lat_c)):
        new_list = model[start:end]
        new_list = np.array(new_list)
        try:
            z =np.vstack((z,new_list))
        except:
            z = [new_list]
            z=np.array(z)
        start+=len(lon_c)
        end+=len(lon_c)

    return z

def model_reciprocal(z,model_indices):
	model = []
	lat_counter = 0
	lon_counter = 0
	for i in range(len(model_indices)):
		lat_i = model_indices[lat_counter,0]
		lon_i = model_indices[lon_counter,1]
		model = np.append(model, z[lat_i,lon_i])
		lat_counter+=1
		lon_counter+=1
	return model
	
def get_tags(obs_refs):

	tag_dict = {'arh':'ANT','ask':'AF','bhd':'OC','bmw':'O','brt':'NA','brw':'NA','cgo':'OC','cmn':'EU','cpt':'AF','cvo':'O','dcc':'ANT','deu':'EU','dig':'EU','glh':'EU','hpb':'EU','irb':'EU','izo':'OC','jfj':'EU', \
				'kmw':'EU','kos':'EU','kvk':'EU','kvv':'EU','lgb':'EU','mlo':'O','mnm':'O','ngl':'EU','nmy':'ANT','nwr':'NA','pal':'EU','pay':'EU','pyr':'AS','rcv':'EU', \
				'rig':'EU','rpb':'O','ryo':'AS','smo':'O','snb':'EU','snl':'SA','spo':'ANT','ssl':'EU','sum':'ARC','syo':'ANT','tkb':'AS','ush':'SA','vdl':'EU', \
				'wes':'EU','yon':'AS','zgt':'EU','zpt':'EU','zrn':'EU'}
	tags = [tag_dict[key] for key in obs_refs]
	tags = np.array(tags)		
		
	return tags
	
def location_process(obs_refs,tags,obs,model):
	loc_types = ['Antarctica','Arctic','Africa','Asia','Europe','Oceania','Oceanic Sites','North America','South America']

	loc_dict_colours


#split data into different arrays by locational type
	for i in range(len(tags)):
		if tags[i] == 'ANT':
			obs_array_0 = np.append(obs_array_0, obs[i])
			model_array_0 = np.append(model_array_0, model[i])
			valid_refs_0 = np.append(valid_refs_0,obs_refs[i])
		if tags[i] == 'ARC':
			obs_array_1 = np.append(obs_array_1, obs[i])
			model_array_1 = np.append(model_array_1, model[i])
			valid_refs_1 = np.append(valid_refs_1,obs_refs[i])
		if tags[i] == 'AF':
			obs_array_2 = np.append(obs_array_2, obs[i])
			model_array_2 = np.append(model_array_2, model[i])
			valid_refs_2 = np.append(valid_refs_2,obs_refs[i])
		if tags[i] == 'AS':
			obs_array_3 = np.append(obs_array_3, obs[i])
			model_array_3 = np.append(model_array_3, model[i])
			valid_refs_3 = np.append(valid_refs_3,obs_refs[i])
		if tags[i] == 'EU':   
			obs_array_4 = np.append(obs_array_4, obs[i])
			model_array_4 = np.append(model_array_4, model[i])
			valid_refs_4 = np.append(valid_refs_4,obs_refs[i])
		if tags[i] == 'OC':
			obs_array_5 = np.append(obs_array_5, obs[i])
			model_array_5 = np.append(model_array_5, model[i])
			valid_refs_5 = np.append(valid_refs_5,obs_refs[i])
		if tags[i] == 'O':
			obs_array_6 = np.append(obs_array_6, obs[i])
			model_array_6 = np.append(model_array_6, model[i])
			valid_refs_6 = np.append(valid_refs_6,obs_refs[i])
		if tags[i] == 'NA':
			obs_array_7 = np.append(obs_array_7, obs[i])
			model_array_7 = np.append(model_array_7, model[i])
			valid_refs_7 = np.append(valid_refs_7,obs_refs[i])
		if tags[i] == 'SA':
			obs_array_8 = np.append(obs_array_8, obs[i])
			model_array_8 = np.append(model_array_8, model[i])
			valid_refs_8 = np.append(valid_refs_8,obs_refs[i])
		
	big_obs = [obs_array_0]+[obs_array_1]+[obs_array_2]+[obs_array_3]+[obs_array_4]+[obs_array_5]+[obs_array_6]+[obs_array_7]+[obs_array_8]#+[obs_array_9]
	big_model = [model_array_0]+[model_array_1]+[model_array_2]+[model_array_3]+[model_array_4]+[model_array_5]+[model_array_6]+[model_array_7]+[model_array_8]#+[model_array_9]
	valid_refs = [valid_refs_0]+[valid_refs_1]+[valid_refs_2]+[valid_refs_3]+[valid_refs_4]+[valid_refs_5]+[valid_refs_6]+[valid_refs_7]+[valid_refs_8]	
	valid_refs = [item for sublist in valid_refs for item in sublist]

	flat_big_obs = [val for subl in big_obs for val in subl]
	flat_big_model = [val for subl in big_model for val in subl]
			
	return big_obs,big_model,valid_refs, flat_big_obs, flat_big_model, loc_types
	
def annual_phase_shift(obs,model):
	obs_test = obs < 9
	obs[obs_test] = obs[obs_test]+12
	model_test= model < 9
	model[model_test] = model[model_test]+12   
	return obs,model
	
def convert_phase_units_single(obs,con_num):			
	for num in range(len(obs)):
		if obs[num] < 0:
			obs[num] = np.pi + (np.pi-np.abs(obs[num]))
		elif obs[num] > 0:
			obs[num] = obs[num]

	convert_factor = con_num/(2*np.pi)
	obs = obs*convert_factor
	
	return obs
    
def convert_phase_units_actual_single(obs,con_num):			
    if obs < 0:
        obs = np.pi+ (np.pi-np.abs(obs))
    elif obs > 0:
		obs = obs

    convert_factor = con_num/(2*np.pi)
    obs = obs*convert_factor
	
    return obs
	
def correct_select_daily_phase_single(obs,lon_c,lat_c,lon_e,lat_e, obs_lon):
    correction = (np.pi*2)/len(lon_c)

#correct obs daily phase output
#for each value, produce a correction number depending on which longitude box the value is in and then correct  
    for num in range(len(obs)):
        lon_i = np.searchsorted(lon_e,obs_lon,side='left')
        lon_i-=1
        if lon_i <len(lon_c)/2:
            lon_i = lon_i+len(lon_c)/2
        elif lon_i >=len(lon_c)/2:
            lon_i = lon_i-len(lon_c)/2
        obs[num] = obs[num]-(lon_i*correction)
        if obs[num] < (-np.pi):
            val = np.abs(obs[num])
            remain = val - np.pi
            val = np.pi - remain
            obs[num] = val
        elif obs[num] > (np.pi):
            val = np.abs(obs[num])
            remain = val - np.pi
            val = -np.pi + remain
            obs[num] = val
    return obs
    
def correct_select_daily_phase_actual_single(obs,lon_c,lat_c,lon_e,lat_e, obs_lon):
    correction = (np.pi*2)/len(lon_c)

#correct obs daily phase output
#for each value, produce a correction number depending on which longitude box the value is in and then correct  
    
    lon_i = np.searchsorted(lon_e,obs_lon,side='left')
    lon_i-=1
    if lon_i <len(lon_c)/2:
        lon_i = lon_i+len(lon_c)/2
    elif lon_i >=len(lon_c)/2:
        lon_i = lon_i-len(lon_c)/2
    obs = obs-(lon_i*correction)
    if obs < (-np.pi):
        val = np.abs(obs)
        remain = val - np.pi
        val = np.pi - remain
        obs = val
    elif obs > (np.pi):
        val = np.abs(obs)
        remain = val - np.pi
        val = -np.pi + remain
        obs = val
    return obs
			