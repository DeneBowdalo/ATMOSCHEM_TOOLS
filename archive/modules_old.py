import numpy as np
import csv
import glob
import datetime as datetime
import matplotlib.pyplot as plt
import logging as log
import scipy.stats as stats
import lomb_phase
import lomb_phase_spec
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq, rfft
import collections
from scipy import signal
import numpy.fft as FFT
import scipy.signal
from math import radians, cos, sin, asin, sqrt 
import math
from mpl_toolkits.basemap import Basemap, shiftgrid, interp, maskoceans
import matplotlib as mpl
from netCDF4 import Dataset
import urllib2
import json
from scipy.interpolate import UnivariateSpline
from scipy.optimize import curve_fit
from incf.countryutils import transformations
from geopy.geocoders import Nominatim
import pandas as pd
from netCDF4 import num2date, date2num

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
	
def date_process(date,time,start_year):
	year=(date//10000)
	month=((date-year*10000)//100)
	day=(date-year*10000-month*100)

	hour=time//100
	min=(time-hour*100)

	doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(int(start_year),1,1,0,0,0) \
              for i in range(len(year))]

	processed_dates=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
	
	return np.array(processed_dates)

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
	print 'top y 1', top_y
	print 'top x 1', top_x
	print 'top indices 1', top_indices

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

	print 'top y 2', top_y
	print 'top x 2', top_x
	print 'top indices 2', top_indices

#remove key peaks not within certain range of top peak
	top_peak_index = np.argmax(top_y)
	top_peak = top_y[top_peak_index]
	limit = (top_peak/100)*90 #peaks must be with 75% of top peak or else removed
	limit = top_peak - limit
	remove_peaks = []
	for counter in range(len(top_y)):
		if top_y[counter] < limit:
			remove_peaks.append(counter)

	top_y = np.delete(top_y,remove_peaks)
	top_x = np.delete(top_x,remove_peaks)
	top_indices = np.delete(top_indices,remove_peaks)

	#print remove_indices

	print 'top y 3', top_y
	print 'top x 3', top_x
	print 'top indices 3', top_indices


#start+=1
# sort cut arrays into ascending order by periods
	re_sorted_x_indices = sorted(range(len(top_x)), key=lambda k: top_x[k])
	sig_x = top_x[re_sorted_x_indices]
	sig_y = top_y[re_sorted_x_indices]
	sig_indices = top_indices[re_sorted_x_indices]

	test = np.logical_or(sig_x < 2,sig_x > 200)
	sig_x = sig_x[test]
	sig_y = sig_y[test]
	sig_indices = sig_indices[test]

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

def obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon): 
#convert obs_lon to same grid as model, if neccessary!
    if obs_lon > lon_e[-1]:
        diff = obs_lon - lon_e[-1]
        obs_lon = lon_e[0] + diff

#check which gridboxes each obs_lat & obs_lon convergance point lie in 
    lat_i = np.searchsorted(lat_e,obs_lat,side='left')
    lon_i = np.searchsorted(lon_e,obs_lon,side='left')
    lat_i = lat_i-1
    lon_i = lon_i-1
	
    return lat_i,lon_i
	
def model_grids(type):
    if type == '4x5':
        #lat lon dims for 4x5 grid  
        lat_c = np.arange(-86.,87.,4)       
        lat_c = np.insert(lat_c,0,-89)
        lat_c = np.append(lat_c,89)

        lon_c = np.arange(-180,176,5)

        lat_e = np.arange(-88.,90,4)
        lat_e = np.insert(lat_e,0,-90.)
        lat_e = np.append(lat_e,90.)

        lon_e = np.arange(-182.5,178,5)

    elif type == '2x2.5':
        #lat lon dims for 2x2.5 grid

        lat_c = np.arange(-88.,89.,2)
        lat_c = np.insert(lat_c,0,-89.5)
        lat_c = np.append(lat_c,89.5)

        lon_c = np.arange(-180,178,2.5)

        lat_e = np.arange(-89.,90,2)
        lat_e = np.insert(lat_e,0,-90.)
        lat_e = np.append(lat_e,90.)

        lon_e = np.arange(-181.25,179,2.5)
    
    elif type == 'GFDL':
        #lat lon dims for GFDL 2x2.5 grid]
        
        lat_e = np.arange(-90,91,2)                                                                                                                                                                                                               
        lon_e = np.arange(-180,181,2.5)
           
        lat_c = np.arange(-89,90,2)                                                                                                                                                                                                               
        lon_c = np.arange(-178.75,179,2.5)
    
    elif type == 'GFDL_nudged':
        #lat lon dims for GFDL 2x2.5 grid]
        
        lat_e = np.arange(-90,91,2)                                                                                                                                                                                                               
        lon_e = np.arange(-180,181,2.5)
           
        lat_c = np.arange(-89,90,2)                                                                                                                                                                                                               
        lon_c = np.arange(-178.75,179,2.5)
    
    return lat_c,lat_e,lon_c,lon_e
    
def take_lomb(time,var,OFAC,SAMP_R,w=False,kp=[]):
#Note - Sampling rate can be lower than average nyquist frequency, as when data is gapped nyquist breaks down. Thus if have something that is non-uniformly spaced, just choose a low samp_r.

    if w == True:
        window = signal.hanning(len(var))
        var_mean = np.mean(var)
        var = var - var_mean
        var = var*window
        amp_corr = 1./(sum(window)/len(window))
    else:
        amp_corr = 1.
    
    #set frequencies (equivalent to FFT freqs)
    N = (len(time)*OFAC)
    T = SAMP_R*N
    df = 1./T
    freqs = np.array([df*n if n<N/2 else df*(n-N) for n in range(N)])
    test = freqs > 0
    freqs = freqs[test]
    
	#take lomb
    fb, mag, ph, fr, fi = lomb_phase.lomb(time,var,freqs)

    periods = 1./freqs

    #CORRECT IMAGINARY COMPONENTS TO FFT EQUIVALENTS
    for i in range(len(fi)):
        if fi[i] < 0:
            fi[i] = fi[i]*-1
        elif fi[i] > 0:
            fi[i] = -fi[i]
    	
    
    #if have key periods calculate amp/phase at specific frequencies and impose into full spectrum where closest point is
    if len(kp) > 0:
        s_periods,s_mag,s_ph,s_fr,s_fi = take_lomb_spec(time,var,False,kp)
        for count in range(len(kp)):
            closest_period = min(range(len(periods)), key=lambda i: abs(periods[i]-kp[count]))
            
            periods[closest_period] = s_periods[count]
            mag[closest_period] = s_mag[count]
            ph[closest_period] = s_ph[count]
            fr[closest_period] = s_fr[count]
            fi[closest_period] = s_fi[count]
    
    if w == True:
        mag = mag * amp_corr
    	
    return periods,mag,ph,fr,fi,amp_corr
    
def take_lomb_spec(time,var,w = False,key_periods=[1./4.,1./3.,1./2.,1.,365.25/4.,365.25/3.,365.25/2.,365.25]):
    if w == True:
        window = signal.hanning(len(var))
        var_mean = np.mean(var)
        var = var - var_mean
        var = var*window
        amp_corr = 1./(sum(window)/len(window))
	
    freqs = 1./np.array(key_periods)
    
	#take lomb
    fb, mag, ph, fr, fi = lomb_phase_spec.lomb(time,var,freqs)
    
    
    if w == True:
        mag = mag * amp_corr
        
    periods = 1./freqs
    
    #CORRECT IMAGINARY COMPONENTS TO FFT EQUIVALENTS
    for i in range(len(fi)):
        if fi[i] < 0:
            fi[i] = fi[i]*-1
        elif fi[i] > 0:
            fi[i] = -fi[i]
    	
    return periods,mag,ph,fr,fi
    

def lomb_ifft(key_periods,time,vals,periods,fr,fi,ofac,amp_corr,w=False):
    #note using windowing reduces length of output time series as need to cut data at start and end of time series. This loss is severe for LSP.
    
    #ONLY WORKS ACCURATELY FOR NON-GAPPED DATA 

    closest_periods = []
    
    for period in key_periods:
        #get indices of key periods, setting real and imag parts to zero that are not relevant.
        #Need to do this as LSP ifft of noise is set equivalent to FFT ifft of noise
        #reconstruction only works for significant periodicities
        closest_period = min(range(len(periods)), key=lambda i: abs(periods[i]-period))
        closest_periods.append(closest_period)

    F = [0]*((len(fr)*2)+2)
    
    #set first real value to average, imag to 0 
    F[0] = complex((np.average(vals)*len(vals)*ofac),0)

    #Get reverse real and imaginary values
    rev_fr=np.copy(fr[::-1])
    rev_fi=np.copy(fi[::-1])

    f_index = 1

    #Fill Fourier Spectrum real and imaginary values
    for i in range(len(fr)):
        F[f_index] = complex(fr[i],fi[i])
        f_index+=1

    F[f_index] = complex(0,0)
    f_index+=1
    
    for i in range(len(fr)):
        F[f_index] = complex(rev_fr[i],-rev_fi[i])
        f_index+=1

    F = np.array(F) 
    
    #Take ifft and just take real values
    ifft_ts = np.fft.ifft(F)
    ifft_ts = ifft_ts.astype('float64')

    #cut reconstructed time series if oversampled
    if ofac > 1:
        actual_len = len(time)
        ifft_ts = ifft_ts[:actual_len]
        if w == True:
            window = signal.hanning(len(vals))
            ifft_ts = ifft_ts/window
            #cut off 40 percent of ends if windowed. Yes loss is severe for LSP
            len_cut = int((len(ifft_ts)/100.)*40.)
            ifft_ts = ifft_ts[len_cut:-len_cut]
            time = time[len_cut:-len_cut]
    
    return time,ifft_ts
    
    
def lomb_ifft_spec(time,vals,periods,mag,ph,samp_r):
    #Reconstruction of time series based on known significant components

    pi2 = np.pi*2.
    
    #get max time length
    max_time = time[-1] - time[0]
    
    #get time array based on max time length and sampling frequency.
    t = []
    t = np.arange(time[0],time[-1]+(samp_r/2.),samp_r)
    
    waveform = np.array([0]*len(t))
        
    #put together waveform with significant component amplitudes and phases.
    for i in range(len(mag)):
        waveform = waveform + (mag[i]*(np.cos((pi2*t/periods[i])-(ph[i]))))
    
    waveform = waveform+np.average(vals) 
        
    return t,waveform

    

def take_fft(time,var,ofac,SAMP_R,w=False,kp=[]):
    var_mean = np.mean(var)
    var = var - var_mean
    
    if w == True:
        window = signal.hanning(len(var))
        var = var*window
        amp_corr = 1./(sum(window)/len(window))
    else:
        amp_corr = 1.
    
    #add zero padding, oversampling of n
    orig_time = np.copy(time)
    n_zeros = (len(time)*ofac) - (len(time))
    zeros_array = [0]*n_zeros
    time = np.append(time,zeros_array)
    var = np.append(var,zeros_array)
    
    #set frequencies (with zero padded oversampling
    n_points = len(time)
    N = np.copy(n_points)
    T = SAMP_R*N
    df = 1./T
    fft_freqs = np.array([df*n if n<N/2 else df*(n-N) for n in range(N)])
    fft_array = fft(var,n_points)
    fft_mag = np.abs(fft_array)
    ratio = n_points/len(var)
    valid = fft_freqs > 0
    fft_periods = 1./fft_freqs
    fft_mag = fft_mag/len(orig_time)
    fft_fr = fft_array.real
    fft_fi = fft_array.imag 
    fft_periods,fft_mag,fft_fr,fft_fi = fft_periods[valid],fft_mag[valid],fft_fr[valid],fft_fi[valid]
    fft_mag = fft_mag*2
    fft_ph = np.arctan2(fft_fi,fft_fr)
    
    #if have key periods calculate amp/phase at specific frequencies and impose into full spectrum where closest point is
    if len(kp) > 0:
        s_periods,s_mag,s_ph,s_fr,s_fi = take_lomb_spec(time,var,False,kp)
        for count in range(len(kp)):
            closest_period = min(range(len(fft_periods)), key=lambda i: abs(fft_periods[i]-kp[count]))
            
            fft_periods[closest_period] = s_periods[count]
            fft_mag[closest_period] = s_mag[count]
            fft_ph[closest_period] = s_ph[count]
            fft_fr[closest_period] = s_fr[count]
            fft_fi[closest_period] = s_fi[count]
    
    if w == True:
        fft_mag = fft_mag * amp_corr   


    return fft_periods,fft_mag,fft_ph,fft_fr,fft_fi,fft_array,amp_corr
    
def fft_ifft(time,vals,fr,fi,fft_array,ofac,amp_corr,w=False):  
   
    #Take ifft and just take real values
    ifft_ts = np.fft.ifft(fft_array)
    ifft_ts = ifft_ts.astype('float64')

    #cut reconstructed time series if oversampled
    if ofac > 1:
        actual_len = len(time)
        ifft_ts = ifft_ts[:actual_len]
        
        if w == True:
            window = signal.hanning(len(vals))
            ifft_ts = ifft_ts/window
            #cut off 2 percent of ends if windowed.
            len_cut = int((len(ifft_ts)/100.)*2.)
            ifft_ts = ifft_ts[len_cut:-len_cut]
            time = time[len_cut:-len_cut]
            ifft_ts = ifft_ts + np.average(vals)
    
    return time,ifft_ts
    
    
def just_smooth(s, periods, n):
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

	#now calculate the rest of the values
	periods_counter = n+1
	for i in s[n+1:]:
		tmp = ( (i - ema[j]) * multiplier ) + ema[j]
		j = j + 1
		ema.append(tmp)
		cut_periods.append(periods[periods_counter])
		periods_counter +=1

	return cut_periods,ema

def log_bin_smooth(limit, n,x_array,y_array):
    #to get exact numbers areas of len n, must add 1 to n. 
    #eg 10 log spaced numbers returns 9 areas to smooth between. Must add 1 to get 10 areas
    n+=1
    
    #generate log spaced numbers
    result = [1]
    if n>1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result)<n:
        next_value = result[-1]*ratio
        if next_value - result[-1] >= 1:
            # safe zone. next_value will be a different integer
            result.append(next_value)
        else:
            # problem! same integer. we need to find next_value by artificially incrementing previous value
            result.append(result[-1]+1)
            # recalculate the ratio so that the remaining values will scale correctly
            ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    # round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
    log_spaced =  np.array(map(lambda x: round(x)-1, result),dtype=np.uint64)
    #amend last value to make it size of input array
    print log_spaced
    log_spaced = log_spaced.tolist()
    log_spaced[-1] = np.uint64(limit)
    
    #start smoothing of arrays
    smoothed_x_array = np.array([])
    smoothed_y_array = np.array([])
    
    for i in range(len(log_spaced)-1):
        try:
            start_i+=1
            end_i+=1
            smoothed_x_array = np.append(smoothed_x_array,np.average(x_array[log_spaced[start_i]:log_spaced[end_i]]))
            smoothed_y_array = np.append(smoothed_y_array,np.average(y_array[log_spaced[start_i]:log_spaced[end_i]]))
        except:
            start_i = 0
            end_i = 1    
            smoothed_x_array = np.append(smoothed_x_array,np.average(x_array[log_spaced[start_i]:log_spaced[end_i]]))
            smoothed_y_array = np.append(smoothed_y_array,np.average(y_array[log_spaced[start_i]:log_spaced[end_i]]))
            
    return smoothed_x_array,smoothed_y_array

def grid_obs_centre_convergance(lat_e,lon_e,obs_lats,obs_lons): 
#convert obs_lon to same grid as model, if neccessary!
    for i in range(len(obs_lons)):
        if obs_lons[i] > lon_e[-1]:
            diff = obs_lons[i] - lon_e[-1]
            obs_lons[i] = lon_e[0] + diff
            print obs_lons[i]
            
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
	
	
def convert_phase_rect(obs,model,con_num):
    convert_factor = (2*np.pi)/con_num
    pi2 = np.pi*2.
	
    for num in range(len(model)):
        model[num] = model[num] * convert_factor
        if model[num] > np.pi:
            model[num] = -np.pi + (model[num]-np.pi)
				
    for num in range(len(obs)):
        obs[num] = obs[num] * convert_factor
        if obs[num] > np.pi:
            obs[num] = -np.pi + (obs[num]-np.pi)
	
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

    tag_dict = {'arh':'ANT','ask':'AF','bhd':'OC','bmw':'O','brt':'NA','brw':'NA','cgo':'OC','cmn':'EU','cpt':'AF','cvo':'O','dcc':'ANT','deu':'EU','dig':'EU','glh':'EU','hpb':'EU','irb':'EU','izo':'O','jfj':'EU', \
                'kmw':'EU','kos':'EU','kps':'EU','kvk':'EU','kvv':'EU','lau':'OC','lgb':'EU','mhd':'EU','mlo':'O','mnm':'O','ngl':'EU','nmy':'ANT','nwr':'NA','pal':'EU','pay':'EU','pyr':'AS','rcv':'EU', \
                'rig':'EU','rpb':'O','ryo':'AS','smo':'O','snb':'EU','snl':'SA','spo':'ANT','ssl':'EU','sum':'ARC','syo':'ANT','thd':'NA','tkb':'AS','ush':'SA','vdl':'EU','bkt':'OC', \
                'wes':'EU','yon':'AS','zgt':'EU','zpt':'EU','zrn':'EU','est':'NA','egb':'NA','ela':'NA','sat':'NA','alg':'NA','kej':'NA','bra':'NA','qak172':'NA','eve419':'NA','vin140':'NA','irl141':'NA', \
                'oxf122':'NA','prk134':'NA','lav410':'NA','dcp114':'NA','par107':'NA','thr422':'NA','cat175':'NA','gth161':'NA','bbe401':'NA','bel116':'NA','mac426':'NA','ckt136':'NA','bvl130':'NA','knz184':'NA', \
                'pin414':'NA','wsp144':'NA','ash135':'NA','cad150':'NA','cha467':'NA','hwf187':'NA','cdr119':'NA','voy413':'NA','cvl151':'NA','gas153':'NA','wst109':'NA','cnt169':'NA','sum156':'NA','egb181':'NA', \
                'bwr139':'NA','lrl117':'NA','mev405':'NA','grb411':'NA','cth110':'NA','yel408':'NA','cow137':'NA','pet427':'NA','pnd165':'NA','sal133':'NA','hox148':'NA','vpi120':'NA','can407':'NA','spd111':'NA', \
                'uvl124':'NA','mkg113':'NA','ped108':'NA','mor409':'NA','alc188':'NA','che185':'NA','kef112':'NA','cnd125':'NA','san189':'NA','psu106':'NA','aca416':'NA','grs420':'NA','shn418':'NA','glr468':'NA', \
                'rom206':'NA','mck131':'NA','cdz171':'NA','mck231':'NA','alh157':'NA','yos404':'NA','bft142':'NA','ana115':'NA','are128':'NA','jot403':'NA','esp127':'NA','abt147':'NA','snd152':'NA','sek430':'NA', \
                'rom406':'NA','stk138':'NA','den417':'NA','how132':'NA','pnf126':'NA','wnc429':'NA','grc474':'NA','con186':'NA','gb0014r':'EU','no0043r':'EU','no0056r':'EU','no0042g':'ARC','gb0037r':'EU', \
                'gb0002r':'EU','fi0009r':'EU','gb0006r':'EU','fi0096g':'EU','es0008r':'EU','gb0052r':'EU','gb0038r':'EU','ie0031r':'EU','gb0039r':'EU','gb0035r':'EU','gb0033r':'EU','gb0031r':'EU','fi0022r':'EU', \
                'gb0043r':'EU','gb0036r':'EU','no0039r':'EU','fi0037r':'EU','ie0001r':'EU','gr0001r':'EU','gb0049r':'EU','gb0051r':'EU','gb0013r':'EU','gb0015r':'EU','no0015r':'EU','gb0050r':'EU','fi0017r':'EU', \
                'es0012r':'EU','gb0048r':'EU','no0052r':'EU','gb0045r':'EU','min':'NA','fre':'NA','pkl':'NA','sna':'NA','gos':'NA','cps':'NA','rtr':'NA','bon':'NA','etl':'NA','fsd':'NA','wsa':'NA', \
                'cz0005r':'EU','fr0008r':'EU','fr0009r':'EU','fr0010r':'EU','fr0013r':'EU','fr0014r':'EU','fr0015r':'EU','fr0016r':'EU','fr0017r':'EU','mk0007r':'EU','nl0009r':'EU','ch0001g':'EU','cz0003r':'EU','nl0091r':'EU', \
                '054703':'NA','090602':'NA','093801':'NA','093901':'NA','080677001':'NA','090019003':'NA','150030011':'O','270370020':'NA','270370423':'NA','380570102':'NA','380570124':'NA','460710001':'NA', \
                '010601':'NA','030501':'NA','031001':'NA','040203':'NA','040207':'NA','040302':'NA','040401':'NA','040501':'NA','040601':'NA','040801':'NA','040901':'NA','041201':'NA','041302':'NA','050103':'NA', \
                '050113':'NA','050116':'NA','050119':'NA','050121':'NA','050126':'NA','050129':'NA','050204':'NA','050308':'NA','050310':'NA','050311':'NA','050404':'NA','050504':'NA','050604':'NA','051501':'NA', \
                '052001':'NA','052201':'NA','052301':'NA','052401':'NA','052601':'NA','052801':'NA','053201':'NA','053301':'NA','053401':'NA','053501':'NA','053601':'NA','053701':'NA','053901':'NA','054102':'NA', \
                '054201':'NA','054801':'NA','054901':'NA','055001':'NA','055101':'NA','055301':'NA','055501':'NA','055601':'NA','055701':'NA','060211':'NA','060410':'NA','060428':'NA','060429':'NA','060513':'NA', \
                '060809':'NA','061004':'NA','061104':'NA','061201':'NA','061603':'NA','061802':'NA','062001':'NA','062501':'NA','062601':'NA','063001':'NA','063301':'NA','063701':'NA','064001':'NA','064101':'NA', \
                '064401':'NA','065001':'NA','065101':'NA','065201':'NA','065401':'NA','065901':'NA','066101':'NA','070118':'NA','070203':'NA','080901':'NA','090120':'NA','090222':'NA','090502':'NA','090701':'NA', \
                '090702':'NA','090801':'NA','091001':'NA','091301':'NA','091401':'NA','091501':'NA','091801':'NA','091901':'NA','092001':'NA','092601':'NA','100110':'NA','100118':'NA','100119':'NA','100125':'NA', \
                '100127':'NA','100128':'NA','100132':'NA','100134':'NA','100135':'NA','101003':'NA','101101':'NA','101202':'NA','101301':'NA','101401':'NA','101501':'NA','101701':'NA','102001':'NA','102102':'NA', \
                '102401':'NA','102801':'NA','103302':'NA','104003':'NA','129102':'NA','010510001':'NA','010731009':'NA','010735003':'NA','010972005':'NA','011170004':'NA','011190002':'NA','020680003':'NA', \
                '040038001':'NA','040058001':'NA','040070010':'NA','040128000':'NA','040131010':'NA','040132001':'NA','040134008':'NA','040134011':'NA','040137003':'NA','040137020':'NA','040137021':'NA', \
                '040137022':'NA','040137024':'NA','040139508':'NA','040139702':'NA','040170119':'NA','040190021':'NA','040191018':'NA','040191020':'NA','040213007':'NA','040217001':'NA','040218001':'NA', \
                '050350005':'NA','051010002':'NA','051130003':'NA','051191002':'NA','060070007':'NA','060090001':'NA','060111002':'NA','060131002':'NA','060194001':'NA','060254003':'NA','060254004':'NA', \
                '060270101':'NA','060290007':'NA','060290011':'NA','060310500':'NA','060390004':'NA','060430003':'NA','060470003':'NA','060650008':'NA','060651016':'NA','060652002':'NA','060670011':'NA', \
                '060675003':'NA','060690002':'NA','060690003':'NA','060710005':'NA','060710012':'NA','060711234':'NA','060719002':'NA','060732007':'NA','060794002':'NA','060798006':'NA','060830008':'NA', \
                '060831013':'NA','060831014':'NA','060831018':'NA','060831021':'NA','060831025':'NA','060833001':'NA','060834003':'NA','060852006':'NA','060893003':'NA','060950005':'NA','061070006':'NA', \
                '061070009':'NA','061110009':'NA','061113001':'NA','061130004':'NA','080013001':'NA','080350004':'NA','080410013':'NA','080590006':'NA','080671004':'NA','080677003':'NA','080690007':'NA', \
                '080690011':'NA','080830101':'NA','090050005':'NA','100010002':'NA','100031007':'NA','100031010':'NA','100051003':'NA','120013011':'NA','120030002':'NA','120050006':'NA','120310077':'NA', \
                '120550003':'NA','120570081':'NA','120573002':'NA','120590004':'NA','120860029':'NA','121010005':'NA','121035002':'NA','121290001':'NA','130210012':'NA','130550001':'NA','130850001':'NA', \
                '131510002':'NA','132130003':'NA','132230003':'NA','132470001':'NA','132611001':'NA','170230001':'NA','170491001':'NA','170650002':'NA','170971007':'NA','171170002':'NA','171570001':'NA', \
                '171971011':'NA','190170011':'NA','191370002':'NA','191630014':'NA','191770006':'NA','201070002':'NA','201730001':'NA','201910002':'NA','201950001':'NA','210150003':'NA','210430500':'NA', \
                '210470006':'NA','210610501':'NA','210910012':'NA','211010014':'NA','211390003':'NA','211850004':'NA','212130004':'NA','212218001':'NA','212270008':'NA','220050004':'NA','220170001':'NA', \
                '220190008':'NA','220190009':'NA','220330013':'NA','220470009':'NA','220470012':'NA','220630002':'NA','220770001':'NA','220890003':'NA','220930002':'NA','220950002':'NA','230052003':'NA', \
                '230090102':'NA','230090103':'NA','230130004':'NA','230173001':'NA','230194008':'NA','230290032':'NA','230310038':'NA','240030014':'NA','240090011':'NA','240130001':'NA','240150003':'NA', \
                '240170010':'NA','240230002':'NA','240251001':'NA','240290002':'NA','240430009':'NA','250010002':'NA','250070001':'NA','250154002':'NA','261130001':'NA','261530001':'NA','270750005':'NA', \
                '271370034':'NA','280110001':'NA','280470008':'NA','290370003':'NA','290390001':'NA','290470003':'NA','290470005':'NA','290470006':'NA','290490001':'NA','291130003':'NA','291370001':'NA', \
                '291570001':'NA','291831002':'NA','291831004':'NA','291860005':'NA','291890005':'NA','291890014':'NA','300298001':'NA','311090016':'NA','320010002':'NA','320030022':'NA','320030023':'NA', \
                '320031019':'NA','320330101':'NA','330074001':'NA','330115001':'NA','340071001':'NA','340110007':'NA','340150002':'NA','340190001':'NA','340230011':'NA','340273001':'NA','340290006':'NA', \
                '340315001':'NA','350010029':'NA','350130008':'NA','350130020':'NA','350130022':'NA','350151005':'NA','350171003':'NA','350431001':'NA','350439004':'NA','350450009':'NA','350450018':'NA', \
                '350451005':'NA','360130011':'NA','360270007':'NA','360310002':'NA','360310003':'NA','360337003':'NA','360410005':'NA','360430005':'NA','360450002':'NA','360530006':'NA','360631006':'NA', \
                '360650004':'NA','360750003':'NA','360790005':'NA','360830004':'NA','360910004':'NA','361010003':'NA','361030004':'NA','361111005':'NA','361173001':'NA','370330001':'NA','370650099':'NA', \
                '370670028':'NA','370671008':'NA','370750001':'NA','370870035':'NA','370870036':'NA','371010002':'NA','371090004':'NA','371170001':'NA','371191009':'NA','371450003':'NA','371570099':'NA', \
                '371590021':'NA','380070002':'NA','380130004':'NA','380250003':'NA','380530002':'NA','380570004':'NA','380650002':'NA','390230001':'NA','390271002':'NA','390410002':'NA','390550004':'NA', \
                '390610010':'NA','390830002':'NA','390950034':'NA','390970007':'NA','391090005':'NA','391351001':'NA','391550009':'NA','391550011':'NA','400019009':'NA','400159008':'NA','400219002':'NA', \
                '400370144':'NA','400430860':'NA','400719010':'NA','400871073':'NA','400979014':'NA','401159004':'NA','401210415':'NA','401430174':'NA','420010002':'NA','420070002':'NA','420070005':'NA', \
                '420270100':'NA','420290100':'NA','420334000':'NA','420430401':'NA','420550001':'NA','420590002':'NA','420630004':'NA','420890002':'NA','420990301':'NA','421174000':'NA','421255001':'NA', \
                '450010001':'NA','450150002':'NA','450190046':'NA','450210002':'NA','450250001':'NA','450290002':'NA','450310003':'NA','450370001':'NA','450730001':'NA','450770002':'NA','450790021':'NA', \
                '450791001':'NA','450830009':'NA','460330132':'NA','470010101':'NA','470090101':'NA','470090102':'NA','470370026':'NA','470651011':'NA','470654003':'NA','470890002':'NA','470930021':'NA', \
                '471210104':'NA','471490101':'NA','471550101':'NA','471571004':'NA','471632002':'NA','471650007':'NA','471650101':'NA','471870106':'NA','471890103':'NA','480290052':'NA','480290059':'NA', \
                '480430101':'NA','481210034':'NA','481830001':'NA','482010029':'NA','482030002':'NA','482090614':'NA','482450101':'NA','483670081':'NA','484230007':'NA','484390075':'NA','484530020':'NA', \
                '490037001':'NA','490370101':'NA','490530130':'NA','500030004':'NA','500070007':'NA','510330001':'NA','510410004':'NA','510610002':'NA','511130003':'NA','511390004':'NA','511630003':'NA', \
                '511970002':'NA','530090013':'NA','530530012':'NA','530531010':'NA','550090026':'NA','550210015':'NA','550290004':'NA','550390006':'NA','550410007':'NA','550590019':'NA','550610002':'NA', \
                '550710007':'NA','550730012':'NA','550850004':'NA','550890008':'NA','550890009':'NA','551110007':'NA','551170006':'NA','551250001':'NA','551270005':'NA','560050123':'NA','560050456':'NA', \
                '560350099':'NA','560350100':'NA','560370200':'NA','560391011':'NA','rishiri':'AS','ogasawara':'O','happo':'AS','hedo':'AS','oki':'AS','ijira':'AS','yusuhara':'AS','alt':'ARC', \
                '040103':'NA','040701':'NA','060104':'NA','060204':'NA','100111':'NA','100202':'NA','040133003':'NA','040134003':'NA','040139997':'NA','040190002':'NA','040191028':'NA','040213003':'NA', \
                '051190007':'NA','060050002':'NA','060070002':'NA','060130002':'NA','060190007':'NA','060195001':'NA','060290008':'NA','060290014':'NA','060290232':'NA','060295001':'NA','060296001':'NA', \
                '060333001':'NA','060370002':'NA','060370016':'NA','060370113':'NA','060371002':'NA','060371103':'NA','060371201':'NA','060371701':'NA','060372005':'NA','060374002':'NA','060375005':'NA', \
                '060376012':'NA','060379033':'NA','060410001':'NA','060450008':'NA','060530002':'NA','060531003':'NA','060550003':'NA','060570005':'NA','060590007':'NA','060591003':'NA','060592022':'NA', \
                '060595001':'NA','060610002':'NA','060610004':'NA','060610006':'NA','060650012':'NA','060655001':'NA','060656001':'NA','060670002':'NA','060670010':'NA','060792006':'NA','060793001':'NA', \
                '060971003':'NA','061111004':'NA','061112002':'NA','061131003':'NA','080050002':'NA','080410016':'NA','080590005':'NA','080691004':'NA','100031013':'NA','120118002':'NA','120210004':'NA', \
                '120330004':'NA','120571035':'NA','120571065':'NA','120830003':'NA','120830004':'NA','120860027':'NA','121012001':'NA','121030004':'NA','121111002':'NA','121130015':'NA','121151005':'NA', \
                '121171002':'NA','121272001':'NA','121275002':'NA','170310076':'NA','170314007':'NA','170314201':'NA','170317002':'NA','170436001':'NA','171193007':'NA','171431001':'NA','171613002':'NA', \
                '171630010':'NA','220330009':'NA','220331001':'NA','220470007':'NA','220511001':'NA','220550007':'NA','220570004':'NA','220730004':'NA','240053001':'NA','250170009':'NA','250250042':'NA', \
                '320030020':'NA','320030043':'NA','320030071':'NA','320030072':'NA','320030073':'NA','320030075':'NA','320030538':'NA','320030601':'NA','320032002':'NA','320310020':'NA','320310025':'NA', \
                '320311005':'NA','320312009':'NA','340170006':'NA','340210005':'NA','350010023':'NA','350010024':'NA','350010027':'NA','350011012':'NA','350011013':'NA','361192004':'NA','371190041':'NA', \
                '371830014':'NA','401091037':'NA','401430137':'NA','420031005':'NA','450790007':'NA','480610006':'NA','480850005':'NA','481130087':'NA','482010026':'NA','482010046':'NA','482010051':'NA', \
                '482010055':'NA','482010062':'NA','482010066':'NA','482011035':'NA','482011039':'NA','482210001':'NA','482311006':'NA','482450009':'NA','482450011':'NA','482450022':'NA','482450102':'NA', \
                '482450628':'NA','482510003':'NA','483390078':'NA','483550025':'NA','483550026':'NA','483611001':'NA','483611100':'NA','483970001':'NA','484391002':'NA','484392003':'NA','484393009':'NA', \
                '490353006':'NA','510595001':'NA','551091002':'NA','800020001':'NA','at0002r':'EU','at0005r':'EU','at0030r':'EU','at0032r':'EU','at0034g':'EU','at0037r':'EU','at0038r':'EU','at0040r':'EU', \
                'at0041r':'EU','at0042r':'EU','at0043r':'EU','at0044r':'EU','at0045r':'EU','at0046r':'EU','at0047r':'EU','at0048r':'EU','be0001r':'EU','be0032r':'EU','be0035r':'EU','bg0053r':'EU', \
                'ch0002r':'EU','ch0003r':'EU','ch0004r':'EU','cy0002r':'EU','cz0001r':'EU','de0001r':'EU','de0002r':'EU','de0003r':'EU','de0007r':'EU','de0008r':'EU','de0009r':'EU','dk0005r':'EU', \
                'dk0031r':'EU','dk0041r':'EU','ee0009r':'EU','ee0011r':'EU','es0007r':'EU','es0009r':'EU','es0010r':'EU','es0011r':'EU','es0013r':'EU','es0014r':'EU','es0016r':'EU','hu0002r':'EU', \
                'it0001r':'EU','lt0015r':'EU','lv0010r':'EU','no0001r':'EU','no0055r':'EU','pl0002r':'EU','pl0003r':'EU','pl0004r':'EU','pl0005r':'EU','pt0004r':'EU','se0005r':'EU','se0011r':'EU', \
                'se0012r':'EU','se0013r':'EU','se0014r':'EU','se0032r':'EU','se0035r':'EU','se0039r':'EU','si0008r':'EU','si0031r':'EU','si0032r':'EU','sk0002r':'EU','sk0004r':'EU','sk0006r':'EU', \
                'sk0007r':'EU','ice':'EU','lyk123':'NA','041101':'NA','052701':'NA','065301':'NA','065701':'NA','091601':'NA','093101':'NA','100314':'NA','060250006':'NA','060731006':'NA','060870003':'NA',\
                '110010025':'NA','120972002':'NA','121152002':'NA','220190002':'NA','350130021':'NA','380171004':'NA','400170101':'NA','482570005':'NA','550030010':'NA','550270007':'NA','550790026':'NA', \
                '800020005':'NA','800020018':'NA','banryu':'AS','ch004':'SA','ch006':'SA','ch007':'SA','350290003':'NA','481211032':'NA','481390016':'NA','si0033r':'EU','nl0010r':'EU','054401':'NA', \
                'abat0enk1':'EU','abat0ill1':'EU','abat0pil1':'EU','abat0zoe2':'EU','abat10002':'EU','abat2f202':'EU','abat2m121':'EU','abat2sp10':'EU','abat2sp20':'EU','abat2vk26':'EU','abat2wo35':'EU', \
                'abat30103':'EU','abat30202':'EU','abat30302':'EU','abat30403':'EU','abat30502':'EU','abat30603':'EU','abat30701':'EU','abat30801':'EU','abat31102':'EU','abat31204':'EU','abat31502':'EU', \
                'abat31701':'EU','abat31904':'EU','abat32101':'EU','abat4s108':'EU','abat4s165':'EU','abat4s418':'EU','abat4s420':'EU','abat52100':'EU','abat53055':'EU','abat56071':'EU','abat60137':'EU', \
                'abat60151':'EU','abat60157':'EU','abat60185':'EU','abat60190':'EU','abat72538':'EU','abat72705':'EU','abat82801':'EU','abbetn016':'EU','abbetn029':'EU','abbetn040':'EU','abbetn046':'EU', \
                'abbetn051':'EU','abbetn052':'EU','abbetn054':'EU','abbetn063':'EU','abbetn073':'EU','abbetn085':'EU','abbetn093':'EU','abbetn100':'EU','abbetn113':'EU','abbetn121':'EU','abbetn132':'EU', \
                'abbetr222':'EU','abbg0055a':'EU','abch0002r':'EU','abch0003r':'EU','abch0016a':'EU','abch0019a':'EU','abch0020a':'EU','abch0024a':'EU','abch0028a':'EU','abch0031a':'EU','abch0033a':'EU', \
                'abch0037a':'EU','abch0040a':'EU','abch0041a':'EU','abch0045a':'EU','abch0051a':'EU','abcy0002r':'EU','abcz0avyn':'EU','abcz0bkuc':'EU','abcz0bmis':'EU','abcz0chvo':'EU','abcz0ckoc':'EU', \
                'abcz0esvr':'EU','abcz0jkmy':'EU','abcz0jkos':'EU','abcz0kprb':'EU','abcz0lsou':'EU','abcz0mjes':'EU','abcz0ppla':'EU','abcz0pprm':'EU','abcz0sonr':'EU','abcz0tbkr':'EU','abcz0tcer':'EU', \
                'abcz0tstd':'EU','abcz0ulom':'EU','abcz0urvh':'EU','abcz0usnz':'EU','abcz0utus':'EU','abcz0uval':'EU','abcz0zsnv':'EU','abdebb053':'EU','abdebb065':'EU','abdebb066':'EU','abdebe027':'EU', \
                'abdebe032':'EU','abdebe056':'EU','abdebe062':'EU','abdebw004':'EU','abdebw030':'EU','abdebw031':'EU','abdebw087':'EU','abdebw103':'EU','abdeby013':'EU','abdeby047':'EU','abdeby049':'EU', \
                'abdeby072':'EU','abdeby081':'EU','abdeby109':'EU','abdehe023':'EU','abdehe024':'EU','abdehe026':'EU','abdehe028':'EU','abdehe042':'EU','abdehe043':'EU','abdehe046':'EU','abdehe051':'EU', \
                'abdehe052':'EU','abdehe060':'EU','abdemv004':'EU','abdemv012':'EU','abdemv017':'EU','abdeni019':'EU','abdeni031':'EU','abdeni051':'EU','abdeni058':'EU','abdeni059':'EU','abdeni060':'EU', \
                'abdeni063':'EU','abdenw064':'EU','abdenw065':'EU','abdenw068':'EU','abdenw074':'EU','abdenw081':'EU','abderp013':'EU','abderp014':'EU','abderp015':'EU','abderp016':'EU','abderp017':'EU', \
                'abderp028':'EU','abdesh001':'EU','abdesh008':'EU','abdesh011':'EU','abdesh017':'EU','abdesl019':'EU','abdesn049':'EU','abdesn051':'EU','abdesn052':'EU','abdesn074':'EU','abdesn076':'EU', \
                'abdesn079':'EU','abdesn080':'EU','abdest068':'EU','abdest069':'EU','abdest089':'EU','abdest098':'EU','abdeth026':'EU','abdeth027':'EU','abdeth040':'EU','abdeth042':'EU','abdeth061':'EU', \
                'abdeub001':'EU','abdeub005':'EU','abdeub028':'EU','abdeub029':'EU','abdeub030':'EU','abdk0031r':'EU','abdk0041a':'EU','abdk0054a':'EU','abee0009r':'EU','abee0011r':'EU','abee0016a':'EU', \
                'abee0020a':'EU','abes0008r':'EU','abes0010r':'EU','abes0011r':'EU','abes0012r':'EU','abes0013r':'EU','abes0014r':'EU','abes0016r':'EU','abes0094a':'EU','abes0316a':'EU','abes0324a':'EU', \
                'abes0330a':'EU','abes0774a':'EU','abes0813a':'EU','abes1173a':'EU','abes1201a':'EU','abes1217a':'EU','abes1222a':'EU','abes1248a':'EU','abes1285a':'EU','abes1311a':'EU','abes1347a':'EU', \
                'abes1359a':'EU','abes1378a':'EU','abes1379a':'EU','abes1400a':'EU','abes1488a':'EU','abes1489a':'EU','abes1491a':'EU','abes1517a':'EU','abes1518a':'EU','abes1531a':'EU','abes1542a':'EU', \
                'abes1543a':'EU','abes1588a':'EU','abes1599a':'EU','abes1616a':'EU','abes1648a':'EU','abes1654a':'EU','abes1660a':'EU','abes1661a':'EU','abes1670a':'EU','abes1671a':'EU','abes1688a':'EU', \
                'abes1689a':'EU','abes1690a':'EU','abes1746a':'EU','abes1753a':'EU','abes1754a':'EU','abes1779a':'EU','abes1793a':'EU','abes1794a':'EU','abfi00208':'EU','abfi00293':'EU','abfi00349':'EU', \
                'abfi00351':'EU','abfi00352':'EU','abfi00356':'EU','abfi00368':'EU','abfi00372':'EU','abfi00428':'EU','abfi00453':'EU','abfi00460':'EU','abfi00564':'EU','abfr02023':'EU','abfr02024':'EU', \
                'abfr03027':'EU','abfr03031':'EU','abfr04038':'EU','abfr04066':'EU','abfr04142':'EU','abfr04158':'EU','abfr04322':'EU','abfr05053':'EU','abfr07020':'EU','abfr07022':'EU','abfr08204':'EU', \
                'abfr08209':'EU','abfr09022':'EU','abfr12020':'EU','abfr12029':'EU','abfr12031':'EU','abfr13011':'EU','abfr14008':'EU','abfr15001':'EU','abfr15012':'EU','abfr16017':'EU','abfr16031':'EU', \
                'abfr18010':'EU','abfr18026':'EU','abfr18037':'EU','abfr18039':'EU','abfr18045':'EU','abfr19001':'EU','abfr20049':'EU','abfr22017':'EU','abfr23124':'EU','abfr24023':'EU','abfr25045':'EU', \
                'abfr26012':'EU','abfr30033':'EU','abfr31008':'EU','abfr32008':'EU','abfr34003':'EU','abfr34043':'EU','abfr34054':'EU','abfr35012':'EU','abfr36005':'EU','abgb0002r':'EU','abgb0006r':'EU', \
                'abgb0013r':'EU','abgb0014r':'EU','abgb0015r':'EU','abgb0031r':'EU','abgb0033r':'EU','abgb0035r':'EU','abgb0036r':'EU','abgb0037r':'EU','abgb0038r':'EU','abgb0039r':'EU','abgb0051a':'EU', \
                'abgr0110r':'EU','abhu0002r':'EU','abhu0040a':'EU','abie0001r':'EU','abie0031r':'EU','abie0090a':'EU','abie0091a':'EU','abie0102a':'EU','abie0111a':'EU','abis0007a':'EU','abit0741a':'EU', \
                'abit0824a':'EU','abit0842a':'EU','abit0952a':'EU','abit0988a':'EU','abit0989a':'EU','abit0992a':'EU','abit1121a':'EU','abit1174a':'EU','abit1179a':'EU','abit1188a':'EU','abit1233a':'EU', \
                'abit1236a':'EU','abit1371a':'EU','abit1373a':'EU','abit1451a':'EU','abit1474a':'EU','abit1519a':'EU','abit1665a':'EU','abit1678a':'EU','abit1685a':'EU','abit1695a':'EU','abli0002a':'EU', \
                'ablt00051':'EU','ablt00053':'EU','ablt00054':'EU','ablu0102a':'EU','ablu0103a':'EU','ablu0104a':'EU','ablu0105a':'EU','abnl00107':'EU','abnl00131':'EU','abnl00227':'EU','abnl00230':'EU', \
                'abnl00235':'EU','abnl00301':'EU','abnl00318':'EU','abnl00437':'EU','abnl00444':'EU','abnl00538':'EU','abnl00620':'EU','abnl00631':'EU','abnl00633':'EU','abnl00641':'EU','abnl00722':'EU', \
                'abnl00738':'EU','abnl00807':'EU','abnl00818':'EU','abnl00918':'EU','abnl00929':'EU','abnl00934':'EU','abno0001r':'EU','abno0015r':'EU','abno0039r':'EU','abno0042r':'EU','abno0043r':'EU', \
                'abno0052r':'EU','abno0055r':'EU','abno0056r':'EU','abpl0002r':'EU','abpl0004r':'EU','abpl0005r':'EU','abpl0014a':'EU','abpl0028a':'EU','abpl0077a':'EU','abpl0094a':'EU','abpl0105a':'EU', \
                'abpl0121a':'EU','abpl0128a':'EU','abpl0150a':'EU','abpl0182a':'EU','abpl0211a':'EU','abpl0243a':'EU','abpl0246a':'EU','abpl0247a':'EU','abpt01044':'EU','abpt02019':'EU','abpt02020':'EU', \
                'abpt03091':'EU','abpt03092':'EU','abpt03096':'EU','abpt04003':'EU','abpt04006':'EU','abro0072a':'EU','abse0005r':'EU','abse0011r':'EU','abse0012r':'EU','abse0013r':'EU','abse0014r':'EU', \
                'abse0026a':'EU','abse0035r':'EU','abse0039r':'EU','abse0054a':'EU','absi0008r':'EU','absi0033a':'EU','absk0004r':'EU','absk0006r':'EU','absk0007r':'EU','absk0041a':'EU','abfi00357':'EU', \
                'abat0ach1':'EU','abat11002':'EU','abat30102':'EU','abat32201':'EU','abat60147':'EU','abat82707':'EU','abdehe002':'EU','abdehe027':'EU','abdehe034':'EU','abdemv001':'EU','abdenw063':'EU', \
                'abdenw075':'EU','abdenw076':'EU','abdesl008':'EU','abdeub002':'EU','abdeub006':'EU','abdeub007':'EU','abdeub012':'EU','abdeub013':'EU','abdeub019':'EU','abdeub022':'EU','abdeub023':'EU', \
                'abdeub024':'EU','abdeub025':'EU','abdeub026':'EU','abdeub027':'EU','abdeub031':'EU','abdeub032':'EU','abdeub033':'EU','abdeub035':'EU','abfi00350':'EU','abnl00724':'EU','abnl00733':'EU', \
                'abnl00913':'EU','abnl00928':'EU','at0004r':'EU','ca0420g':'ARC','de0004r':'EU','de0011r':'EU','de0012r':'EU','de0017r':'EU','de0026r':'EU','de0035r':'EU','es0001r':'EU','es0003r':'EU',  \
                'es0004r':'EU','es0005r':'EU','fi0004r':'EU','it0004r':'EU','no0030r':'EU','no0041r':'EU','no0044r':'EU','no0045r':'EU','no0048r':'EU','no0488r':'EU','se0002r':'EU','cha':'NA','wel149':'NA', \
                '051901':'NA','052101':'NA','052901':'NA','053001':'NA','053101':'NA','053801':'NA','060302':'NA','061005':'NA','062201':'NA','062401':'NA','062701':'NA','063601':'NA','064201':'NA','090901':'NA', \
                '050970001':'NA','060530006':'NA','060675002':'NA','060770009':'NA','060831012':'NA','060831015':'NA','060831016':'NA','060831019':'NA','061010002':'NA','061110004':'NA','061110005':'NA','120813002':'NA', \
                '120860030':'NA','121091003':'NA','191131015':'NA','191632011':'NA','210610500':'NA','220110002':'NA','220190007':'NA','220470002':'NA','230090101':'NA','340010005':'NA','450870001':'NA', \
                'edm':'EU','zep':'EU','de0031r':'EU','054301':'NA','380910001':'NA','abdesh013':'EU','abdesh014':'EU','abdeub021':'EU','abnl00232':'EU','ash235':'NA','dev412':'NA','060773003':'NA','060792004':'NA', \
                '220430001':'NA','abat60184':'EU','abcz0uvse':'EU','abdesn057':'EU','gb0044r':'EU','ncs415':'NA','030801':'NA','530570013':'NA','abat30992':'EU','abat31496':'EU','abcz0mbup':'EU','abdehe039':'EU', \
                'abdest070':'EU','abgb0044r':'EU','abgb0617a':'EU','mt0001r':'EU','060832012':'NA','450110001':'NA','abdeub034':'EU','abes1200a':'EU','abes1436a':'EU','abgb0045r':'EU','abpl0024a':'EU', \
                'abdeub038':'EU','abdeub039':'EU','abgb0043r':'EU','ablv00010':'EU','abno0041r':'EU','abno0045r':'EU','abpl0026a':'EU','de0039r':'EU','vii423':'O','780200001':'O','abdenw093':'EU', \
                'abdeub017':'EU','abdeub040':'EU','abes1216a':'EU','abfr04324':'EU','abfr04328':'EU','abfr16201':'EU','abfr19009':'EU','abfr22010':'EU','abit0774a':'EU','abpl0030a':'EU','de0045r':'EU', \
                'de0047r':'EU','030701':'NA','150010006':'O','tappi':'AS','abes1538a':'EU','absk0008a':'EU','absk0009a':'EU','aht':'EU','oul':'EU','uto':'EU','vir':'EU','120570110':'NA','abgb0745a':'EU', \
                'abes1605a':'EU','abes1607a':'EU','abes1614a':'EU','abgb0754a':'EU','030901':'NA','055201':'NA','460711001':'NA','abat30407':'EU','abit1522a':'EU','abro0060a':'EU','nl0007r':'EU','abbg0056a':'EU', \
                'abes1662a':'EU','abes1802a':'EU','abes1805a':'EU','abes1806a':'EU','abes1808a':'EU','abes1810a':'EU','abes1811a':'EU','abes1813a':'EU','abfr10041':'EU','abfr23175':'EU','abfr24030':'EU','abfr26019':'EU', \
                'abfr34034':'EU','abit1464a':'EU','abit1681a':'EU','abit1736a':'EU','ablt00052':'EU','abpt02021':'EU','abro0082a':'EU','060303':'NA','066201':'NA','091101':'NA','060798005':'NA','481391044':'NA', \
                '483091037':'NA','abes1827a':'EU','abes1851a':'EU','abfr07038':'EU','abfr12046':'EU','abfr25049':'EU','abfr38011':'O','abfr38012':'O','abfr41007':'EU','abit0659a':'EU','abit1061a':'EU','abit1149a':'EU', \
                'abit1214a':'EU','abit1288a':'EU','abit1340a':'EU','abit1553a':'EU','abit1619a':'EU','abit1848a':'EU','abit1863a':'EU','abit1865a':'EU','abmt00007':'EU','abpt01047':'EU','abse0032r':'EU','abse0066a':'EU', \
                'stn012':'AS','stn014':'AS','stn018':'ARC','stn021':'NA','stn024':'ARC','stn029':'O','stn043':'O','stn053':'EU','stn076':'NA','stn077':'NA','stn089':'ARC','stn099':'EU','stn101':'ANT','stn107':'NA','stn109':'O', \
                'stn156':'EU','stn175':'AF','stn191':'O','stn233':'ANT','stn308':'EU','stn315':'ARC','stn316':'EU','stn318':'EU','stn323':'ANT','stn330':'AS','stn338':'NA','stn348':'EU','stn394':'OC','stn435':'SA','stn436':'O', \
                'stn437':'OC','stn443':'OC','stn456':'NA','stn457':'NA','stn458':'NA','ascen':'O','costarica':'SA','hanoi':'AS','hilo':'O','java':'OC','kuala':'OC','nairobi':'AF','natal':'SA','paramaribo':'SA','reunion':'O', \
                'samoa':'O','ktb':'EU','mhn':'EU','nia':'EU','roq':'EU','spm':'EU','be0011r':'EU','be0013r':'EU','de0043g':'EU','de0044r':'EU','dk0008r':'EU','dk0009r':'EU','dk0012r':'EU','es0006r':'EU','es0017r':'EU', \
                'gb0053r':'EU','lv0016r':'EU','md0013r':'EU','nl0002r':'EU','nl0011r':'EU','no0002r':'EU','no0008r':'EU','no0047r':'EU','no0057r':'EU','no0091r':'EU','no0099r':'EU','no0796r':'EU','no0797r':'EU', \
                'no1010r':'EU','ro0008r':'EU','rs0005r':'EU','ru0001r':'EU','ru0013r':'EU','ru0014r':'EU','ru0016r':'EU','se0008r':'EU','sk0005r':'EU','080801':'NA','090606':'NA','090805':'NA','090806':'NA', \
                '092201':'NA','092401':'NA','092501':'NA','092901':'NA','093401':'NA','093501':'NA','093701':'NA','094201':'NA','094202':'NA','094401':'NA','100315':'NA','101803':'NA','104301':'NA','105301':'NA', \
                '105601':'NA','105604':'NA','106800':'NA','129103':'NA','129202':'NA','010710020':'NA','010830004':'NA','021221004':'NA','040139993':'NA','060012005':'NA','060231005':'NA','060610007':'NA','060773002':'NA', \
                '060791004':'NA','060794001':'NA','060831011':'NA','060831017':'NA','060831026':'NA','060831027':'NA','060831030':'NA','060834004':'NA','060835001':'NA','060991005':'NA','061110006':'NA','090131001':'NA', \
                '160210003':'NA','171050002':'NA','180050007':'NA','180330002':'NA','180510010':'NA','180730003':'NA','181090004':'NA','181290001':'NA','181290002':'NA','181410012':'NA','181470002':'NA','181470006':'NA', \
                '181470008':'NA','181550001':'NA','181830003':'NA','211771004':'NA','220330008':'NA','220470006':'NA','260050001':'NA','260050002':'NA','260210013':'NA','260430901':'NA','260430902':'NA','261050005':'NA', \
                '261110941':'NA','261390006':'NA','270177416':'NA','270495302':'NA','270757630':'NA','270757631':'NA','270757632':'NA','270757634':'NA','270757635':'NA','271710007':'NA','280030004':'NA','280190001':'NA', \
                '280450001':'NA','280930001':'NA','281070001':'NA','291290001':'NA','291830010':'NA','292110001':'NA','300710010':'NA','300750001':'NA','300830001':'NA','300870001':'NA','301110086':'NA','330150013':'NA', \
                '350015010':'NA','350281002':'NA','350451233':'NA','370370004':'NA','371050002':'NA','380130001':'NA','380130002':'NA','380570103':'NA','380570104':'NA','381010114':'NA','390010013':'NA','391030004':'NA', \
                '391291001':'NA','391450014':'NA','400719003':'NA','400770441':'NA','401010160':'NA','401010167':'NA','401110152':'NA','401110153':'NA','401359015':'NA','410591003':'NA','420010001':'NA','420150011':'NA', \
                '440030002':'NA','450750003':'NA','450791006':'NA','461094003':'NA','461270001':'NA','461270002':'NA','470110004':'NA','470310004':'NA','470430009':'NA','470630003':'NA','470750002':'NA','470750003':'NA', \
                '470850020':'NA','471050106':'NA','471190106':'NA','471250009':'NA','471251010':'NA','471310004':'NA','471410004':'NA','471451020':'NA','471572005':'NA','481390017':'NA','483491051':'NA','490530007':'NA', \
                '511650002':'NA','530150014':'NA','530330023':'NA','540250001':'NA','540990003':'NA','550210016':'NA','560090819':'NA','cheju':'AS','jinyunshan':'AS','kanghwa':'AS','abat0hbg1':'EU','abat2m226':'EU', \
                'abat2sv24':'EU','abat30402':'EU','abat30602':'EU','abat30901':'EU','abat31703':'EU','abat31902':'EU','abat31903':'EU','abat31905':'EU','abat31906':'EU','abat32604':'EU','abat4s199':'EU','abat4s407':'EU', \
                'abat60103':'EU','abat60104':'EU','abat60114':'EU','abat60115':'EU','abat60119':'EU','abat60144':'EU','abat60154':'EU','abat60180':'EU','abat60197':'EU','abba0001g':'EU','abbelhr01':'EU','abbelld02':'EU', \
                'abbetn015':'EU','abbetn027':'EU','abbetn060':'EU','abbg0019a':'EU','abbg0026a':'EU','abbg0038a':'EU','abbg0057a':'EU','abbg0058a':'EU','abbg0071a':'EU','abcz0htrm':'EU','abcz0jtre':'EU','abcz0kchm':'EU', \
                'abcz0knal':'EU','abcz0lbli':'EU','abcz0lclm':'EU','abcz0lfru':'EU','abcz0ljnm':'EU','abcz0molj':'EU','abcz0molo':'EU','abcz0ppls':'EU','abcz0sdub':'EU','abcz0skls':'EU','abcz0skry':'EU','abcz0tbom':'EU', \
                'abcz0thar':'EU','abcz0tora':'EU','abcz0torv':'EU','abcz0tovk':'EU','abcz0tver':'EU','abcz0uchm':'EU','abcz0ukru':'EU','abcz0umed':'EU','abcz0utpm':'EU','abdebb010':'EU','abdebb016':'EU','abdebb037':'EU', \
                'abdebb051':'EU','abdeby030':'EU','abdeby067':'EU','abdeby069':'EU','abdeby0691':'EU','abdeby092':'EU','abdeby124':'EU','abdehe048':'EU','abdehe050':'EU','abdehh063':'EU','abdemv024':'EU','abdeni014':'EU', \
                'abdeni030':'EU','abdeni0311':'EU','abdeni077':'EU','abdenw001':'EU','abdenw005':'EU','abdenw032':'EU','abdenw033':'EU','abdest104':'EU','abdeth017':'EU','abdeth0261':'EU','abdeth037':'EU','abdeub0011':'EU', \
                'abdeub0021':'EU','abdeub0051':'EU','abdeub0061':'EU','abdeub0071':'EU','abdeub008':'EU','abdeub0121':'EU','abdeub0131':'EU','abdeub014':'EU','abdeub0171':'EU','abdeub018':'EU','abdeub036':'EU','abdeub037':'EU', \
                'abdeub041':'EU','abdeub042':'EU','abdk0012r':'EU','abdk0031r1':'EU','abdk0048a':'EU','abee0009r1':'EU','abes0001r':'EU','abes0005r':'EU','abes0006r':'EU','abes0017r':'EU','abes0296a':'EU','abes0297a':'EU', \
                'abes1716a':'EU','abes1774a':'EU','abes1778a':'EU','abes1831a':'EU','abes1853a':'EU','abes1854a':'EU','abes1878a':'EU','abes1882a':'EU','abes1898a':'EU','abes1913a':'EU','abes1917a':'EU','abes1923a':'EU', \
                'abes1987a':'EU','abes1990a':'EU','abes1996a':'EU','abes2014a':'EU','abes2018a':'EU','abfi00424':'EU','abfi00582':'EU','abfr01004':'EU','abfr01014':'EU','abfr02026':'EU','abfr05087':'EU','abfr06018':'EU', \
                'abfr06134':'EU','abfr08401':'EU','abfr09302':'EU','abfr10007':'EU','abfr10014':'EU','abfr10029':'EU','abfr12028':'EU','abfr12052':'EU','abfr12053':'EU','abfr14042':'EU','abfr14051':'EU','abfr15002':'EU', \
                'abfr15018':'EU','abfr16054':'EU','abfr16056':'EU','abfr16058':'EU','abfr16302':'EU','abfr17016':'EU','abfr18049':'EU','abfr20039':'EU','abfr20044':'EU','abfr20068':'EU','abfr21004':'EU','abfr21031':'EU', \
                'abfr22004':'EU','abfr22016':'EU','abfr22022':'EU','abfr23096':'EU','abfr23097':'EU','abfr23156':'EU','abfr23177':'EU','abfr24033':'EU','abfr24039':'EU','abfr27001':'EU','abfr27101':'EU','abfr30029':'EU', \
                'abfr30030':'EU','abfr30037':'EU','abfr33302':'EU','abfr35001':'EU','abfr350011':'EU','abfr36019':'EU','abfr39008':'O','abfr41005':'EU','abfr41024':'EU','abgb0737a':'EU','abgb0838a':'EU','abgb0957a':'EU', \
                'abgb0998a':'EU','abgb0999a':'EU','abgb1017a':'EU','abgr0033a':'EU','abgr0120a':'EU','abie0113a':'EU','abie0115a':'EU','abie0120a':'EU','abie0122a':'EU','abie0124a':'EU','abie0130a':'EU','abie0134a':'EU', \
                'abie0138a':'EU','abie0139a':'EU','abie0144a':'EU','abie0146a':'EU','abie0147a':'EU','abie0148a':'EU','abis0016a':'EU','abis0020a':'EU','abit0267a':'EU','abit0440a':'EU','abit0558a':'EU','abit0692a':'EU', \
                'abit0709a':'EU','abit0838a':'EU','abit0864a':'EU','abit0870a':'EU','abit0982a':'EU','abit0988a1':'EU','abit1085a':'EU','abit1184a':'EU','abit1210a':'EU','abit1212a':'EU','abit1220a':'EU','abit1222a':'EU', \
                'abit1225a':'EU','abit1292a':'EU','abit1294a':'EU','abit1307a':'EU','abit1312a':'EU','abit1418a':'EU','abit1459a':'EU','abit1474a1':'EU','abit1542a':'EU','abit1543a':'EU','abit1548a':'EU','abit1582a':'EU', \
                'abit1595a':'EU','abit1596a':'EU','abit1646a':'EU','abit1663a':'EU','abit1666a':'EU','abit1725a':'EU','abit1729a':'EU','abit1773a':'EU','abit1806a':'EU','abit1823a':'EU','abit1831a':'EU','abit1832a':'EU', \
                'abit1861a':'EU','abit1868a':'EU','abit1870a':'EU','abit1898a':'EU','abit1902a':'EU','abit1904a':'EU','abit1911a':'EU','abit1912a':'EU','abit1914a':'EU','abit1915a':'EU','abit1917a':'EU','abit1919a':'EU', \
                'abit1921a':'EU','abit1924a':'EU','abit1925a':'EU','abit1942a':'EU','abit1944a':'EU','abit1948a':'EU','abit1960a':'EU','abit1998a':'EU','abit2011a':'EU','abit2013a':'EU','abit2014a':'EU','abit2023a':'EU', \
                'abit2027a':'EU','abit2032a':'EU','abit2054a':'EU','abit2060a':'EU','abit2063a':'EU','abit2069a':'EU','abit2071a':'EU','abit2072a':'EU','abit2074a':'EU','abit2076a':'EU','abli0001a':'EU','ablt00021':'EU', \
                'ablu0099a':'EU','ablv000o1':'EU','ablv00ng1':'EU','abnl001071':'EU','abnl001311':'EU','abnl002271':'EU','abnl002351':'EU','abnl00245':'EU','abnl00246':'EU','abnl00247':'EU','abnl003011':'EU','abnl003181':'EU', \
                'abnl004371':'EU','abnl005381':'EU','abnl006201':'EU','abnl006311':'EU','abnl006331':'EU','abnl00644':'EU','abnl00703':'EU','abnl007221':'EU','abnl007241':'EU','abnl008071':'EU','abnl009181':'EU','abnl009281':'EU', \
                'abnl009341':'EU','abno0063a':'EU','abpl0027a':'EU','abpl0043a':'EU','abpl0062a':'EU','abpl0068a':'EU','abpl0115a':'EU','abpl0157a':'EU','abpl0157a1':'EU','abpl0191a':'EU','abpl0198a':'EU','abpl0212a':'EU', \
                'abpl0222a':'EU','abpl0273a':'EU','abpl0276a':'EU','abpl0313a':'EU','abpl0314a':'EU','abpl0317a':'EU','abpl0321a':'EU','abpl0349a':'EU','abpl0437a':'EU','abpl0450a':'EU','abpl0468a':'EU','abpl0469a':'EU', \
                'abpl0470a':'EU','abpl0473a':'EU','abpl0503a':'EU','abpl0505a':'EU','abpl0506a':'EU','abpl0525a':'EU','abpl0552a':'EU','abpl0560a':'EU','abpl0561a':'EU','abpl0568a':'EU','abpl0573a':'EU','abpl0581a':'EU', \
                'abpt01051':'EU','abpt01054':'EU','abpt02022':'EU','abpt03099':'EU','abpt03101':'EU','abpt03102':'EU','abpt04002':'EU','abpt05012':'EU','abpt07001':'O','abro0008r':'EU','abro0077a':'EU','abro0087a':'EU', \
                'abro0126a':'EU','abro0138a':'EU','abro0153a':'EU','abro0191a':'EU','abro0198a':'EU','abro0210a':'EU','abrs0005r':'EU','abrs0031a':'EU','absi0046a':'EU','absk0006a':'EU','absk0007a':'EU','absk0010a':'EU', \
                'absk0011a':'EU','absk0013a':'EU','absk0028a':'EU','absk0030a':'EU','absk0031a':'EU','absk0051a':'EU','abdebw0311':'EU','abdest0701':'EU','abes1662a1':'EU','abit1019a':'EU','abpl0397a':'EU','abro0212a':'EU', \
                'cdl':'NA','chm':'NA','esp':'NA','llb':'NA','bur':'EU','ivn':'EU','jcz':'EU','kam':'EU','leb':'EU','log':'EU','plm':'SA','plv':'EU','zsn':'EU','abes1540a':'EU','abhu0017a':'EU','abhu0018a':'EU','abhu0019a':'EU',
                'ablv000101':'EU','ablv00016':'EU','abno0002r':'EU','abno0008r':'EU','abpl0067a':'EU','abpl0081a':'EU','abpl0084a':'EU','abpl0088a':'EU','abpl0090a':'EU','abpl0103a':'EU','abpl0116a':'EU','abpl0131a':'EU',
                'abpl0153a':'EU','abpl0154a':'EU','abpl0155a':'EU','abpl0158a':'EU','abpl0160a':'EU','abpl0161a':'EU','abpl0163a':'EU','abpl0185a':'EU','abpl0201a':'EU','abpl0208a':'EU','abpl0229a':'EU','abpl0230a':'EU',
                'abpl0232a':'EU','abpl0233a':'EU','abpl0255a':'EU','abpl0259a':'EU','abpl0277a':'EU','abpl0281a':'EU','abpl0282a':'EU','abpl0287a':'EU','abpl0323a':'EU','abpl0324a':'EU','abpl0332a':'EU','abpl0338a':'EU',
                'abpl0340a':'EU','abpl0341a':'EU','abpl0342a':'EU','abpl0344a':'EU','abpl0345a':'EU','abpl0347a':'EU','abpl0351a':'EU','abpl0363a':'EU','abpl0366a':'EU','abpl0370a':'EU','abpl0376a':'EU','abpl0378a':'EU',
                'abpl0379a':'EU','abpl0380a':'EU','abpl0384a':'EU','abpl0389a':'EU','abpl0390a':'EU','abpl0391a':'EU','abpl0392a':'EU','abpl0393a':'EU','abpl0394a':'EU','abpl0418a':'EU','abpl0427a':'EU','abpl0428a':'EU',
                'abpl0430a':'EU','abpl0433a':'EU','abpl0435a':'EU','abpl0439a':'EU','abpl0440a':'EU','abpl0442a':'EU','abpl0451a':'EU','abpl0452a':'EU','abpl0453a':'EU','abpl0454a':'EU','abpl0456a':'EU','abpl0459a':'EU',
                'abpl0460a':'EU','abpl0462a':'EU','abpl0464a':'EU','abpl0471a':'EU','abpl0511a':'EU','abro0001a':'EU','abro0006a':'EU','abro0034a':'EU','abro0038a':'EU','abrs0010a':'EU','abrs0011a':'EU','abrs0013a':'EU',
                'abrs0019a':'EU','abrs0034a':'EU','abrs0040a':'EU','abrs0046a':'EU','abse0002r':'EU','abse0008r':'EU','abse0010a':'EU','abse0018a':'EU','abse0061a':'EU','ams':'O','asc':'O','azr':'O','bal':'EU','bme':'O',
                'bsc':'EU','cba':'NA','cfa':'OC','chr':'O','cmo':'NA','cri':'AS','crz':'O','cya':'ANT','eic':'O','gmi':'O','goz':'EU','hba':'ANT','hun':'EU','itn':'NA','kum':'O','kzd':'AS','lef':'NA','lmp':'EU','maa':'ANT',
                'mbc':'ARC','mid':'O','mqa':'O','nmb':'AF','poc900n':'O','poc905n':'O','poc905s':'O','poc910n':'O','poc910s':'O','poc915n':'EU','poc915s':'EU','poc920n':'O','poc920s':'O','poc925n':'O','poc925s':'O',
                'poc930n':'O','poc930s':'O','poc935n':'O','poc935s':'O','psa':'ANT','pta':'NA','scs903n':'AS','scs906n':'AS','scs909n':'AS','scs912n':'AS','scs915n':'AS','scs918n':'AS','scs921n':'AS','sey':'O','sgp':'NA',
                'shm':'O','sis':'O','stm':'O','tap':'AS','tdf':'SA','uum':'AS','wis':'EU','abrs0044a':'EU','090605':'NA','khanchanaburi':'AS','nakhonratchasima':'AS','ochiishi':'AS','sado-seki':'AS','abes1774a1':'EU',
                'abes1793a1':'EU','abat10007':'EU','abbetn041':'EU','abdeby026':'EU','abdeby093':'EU','abdeni0310':'EU','abdest017':'EU','abdest018':'EU','abdest106':'EU','abdeub0010':'EU','abdeub0020':'EU','abdeub0050':'EU',
                'abdeub0060':'EU','abdeub0070':'EU','abdeub0120':'EU','abdeub0130':'EU','abdeub0170':'EU','abes2001a':'EU','abfi00354':'EU','abfi00426':'EU','abfr03085':'EU','abfr03088':'EU','abfr05088':'EU','abfr06133':'EU',
                'abfr08618':'EU','abfr10004':'EU','abfr10132':'EU','abfr1227a':'EU','abfr13001':'EU','abfr19008':'EU','abfr21050':'EU','abfr23102':'EU','abfr23230':'EU','abfr24001':'EU','abfr29437':'EU','abfr34026':'EU',
                'abfr34038':'EU','abgb0039r0':'EU','abgb0041r':'EU','abgb0048r':'EU','abgb0794a':'EU','abgb0881a':'EU','abie0092a':'EU','abie0109a':'EU','abit0988a0':'EU','abit1121a0':'EU','abit1726a':'EU','abit1758a':'EU',
                'abit1759a':'EU','abit1982a':'EU','abit2022a':'EU','ablt000510':'EU','ablt000530':'EU','ablt000540':'EU','ablt0052a':'EU','abmt00001':'EU','abnl001070':'EU','abnl001310':'EU','abnl002270':'EU','abnl002350':'EU',
                'abnl003010':'EU','abnl003180':'EU','abnl004370':'EU','abnl005380':'EU','abnl006200':'EU','abnl006310':'EU','abnl006330':'EU','abnl007220':'EU','abnl007240':'EU','abnl008070':'EU','abnl009180':'EU','abnl009280':'EU',
                'abnl009340':'EU','abno0056r0':'EU','abpl0383a':'EU','abpt040020':'EU','absi0053a':'EU','absk0010r':'EU','010010003':'NA','010270001':'NA','010499991':'NA','010610001':'NA','010790002':'NA','010970025':'NA',
                '011011001':'NA','040190013':'NA','050199991':'NA','060110002':'NA','060150002':'NA','060150003':'NA','060172003':'NA','060192009':'NA','060379034':'NA','060390500':'NA','060410002':'NA','060651004':'NA',
                '060651010':'NA','060651999':'NA','060731201':'NA','061070005':'NA','061070008':'NA','090159991':'NA','120310106':'NA','120350004':'NA','120730002':'NA','120779991':'NA','120990006':'NA','131910002':'NA',
                '132319991':'NA','170190007':'NA','170191001':'NA','170331001':'NA','170859991':'NA','171170001':'NA','171199991':'NA','171332001':'NA','180030001':'NA','180839991':'NA','181630013':'NA','181699991':'NA',
                '191210003':'NA','191530024':'NA','191770004':'NA','201619991':'NA','201730018':'NA','210350004':'NA','211759991':'NA','211770005':'NA','212210001':'NA','212219991':'NA','212299991':'NA','220190006':'NA',
                '230039991':'NA','230090003':'NA','230194007':'NA','230199991':'NA','240190004':'NA','240199991':'NA','240210034':'NA','240451004':'NA','260492001':'NA','260630905':'NA','260770906':'NA','261330901':'NA',
                '261579991':'NA','261611001':'NA','261619991':'NA','261659991':'NA','270031001':'NA','270710101':'NA','271376317':'NA','271636015':'NA','280370001':'NA','280450719':'NA','280890002':'NA','281250001':'NA',
                '281619991':'NA','310550032':'NA','311079991':'NA','320038000':'NA','330099991':'NA','340090003':'NA','360291005':'NA','360310005':'NA','360319991':'NA','361099991':'NA','370230004':'NA','370319991':'NA',
                '371139991':'NA','371230099':'NA','371239991':'NA','380570101':'NA','380650101':'NA','390030002':'NA','390479991':'NA','391219991':'NA','391532004':'NA','400190296':'NA','400190297':'NA','400310649':'NA',
                '400330680':'NA','400670671':'NA','400770440':'NA','400819005':'NA','400819024':'NA','400850300':'NA','400892001':'NA','401430177':'NA','410290008':'NA','420019991':'NA','420279991':'NA','420479991':'NA',
                '420859991':'NA','421119991':'NA','450030004':'NA','450070005':'NA','450230002':'NA','450290001':'NA','450451003':'NA','450890001':'NA','460110003':'NA','461270003':'NA','470250001':'NA','470259991':'NA',
                '470419991':'NA','470550001':'NA','471190016':'NA','480710900':'NA','480710902':'NA','480710903':'NA','483739991':'NA','484570101':'NA','510150004':'NA','510719991':'NA','511479991':'NA','511870002':'NA',
                '530530027':'NA','540219991':'NA','540939991':'NA','550170001':'NA','550210005':'NA','550210008':'NA','550210013':'NA','550270001':'NA','550330003':'NA','550370001':'NA','550430003':'NA','551199991':'NA',
                '551230008':'NA','560370898':'NA','are228':'NA','dcp214':'NA','lcw121':'NA','onl102':'NA','pbf129':'NA','sum256':'NA','wfm105':'NA','wpa103':'NA','wpb104':'NA','020201':'NA','054001':'NA','062101':'NA',
                '063401':'NA','063901':'NA','064301':'NA','066001':'NA','129401':'ARC','129501':'NA','at0003r':'EU','ch0031r':'EU','de0006r':'EU','de0013r':'EU','de0018r':'EU','de0038r':'EU','de0042r':'EU','de0046r':'EU',
                'dk0010g':'ARC','fr0011r':'EU','fr0018r':'EU','fr0032r':'EU','gb0041r':'EU','gr0002r':'EU','no0762r':'EU','ru0018r':'EU','se0003r':'EU','se0033r':'EU','se0034r':'EU','se0094r':'EU','ang':'O','cas':'EU',
                'dak':'EU','dbl':'EU','dmv':'AS','lon':'NA','mbi':'ANT','mcm':'ANT','mvh':'EU','shp':'EU','sja':'SA','wkt':'NA','fra':'NA','mtm':'NA','abdeby0690':'EU','abdeth0260':'EU','abdk0031r0':'EU','abee0009r0':'EU',
                'abit1474a0':'EU','abpl0157a0':'EU','180730002':'NA','180770001':'NA','180830004':'NA','181270015':'NA','181270016':'NA','181291002':'NA','181530001':'NA','181631001':'NA','181631002':'NA','181671012':'NA',
                '210150006':'NA','210150007':'NA','212230004':'NA','261470904':'NA','261470905':'NA','550110004':'NA','550110007':'NA','550350010':'NA','550591001':'NA','551230003':'NA','abdebw0310':'EU','abdest0700':'EU',
                'abes1662a0':'EU','120310086':'NA','abes1774a0':'EU','371630003':'NA','371630004':'NA','abfr27006':'EU','ablv000100':'EU','abfr08025':'EU','abfr20204':'EU','abfr36010':'EU','230230006':'NA','230290019':'NA',
                '230310040':'NA','240338003':'NA','270953051':'NA','290190011':'NA','370510008':'NA','371290002':'NA'}
                
    for i in range(len(obs_refs)):
        if '_' in obs_refs[i]:
            split_res = obs_refs[i].split('_')
            obs_refs[i] = split_res[0]
            
    inv = []
    for i in obs_refs:
        if i not in tag_dict:
            inv.append(i)

    print ("\n%s\n"%("':'NA','".join(i for i in inv)))
        
    tags = [tag_dict[key] for key in obs_refs]
    tags = np.array(tags)		

    return tags
	
def location_process(obs_refs,tags,obs,model):
	loc_types = ['Antarctica','Arctic','Africa','Asia','Europe','Oceania','Oceanic Sites','North America','South America']

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
	obs_test = obs < 6
	obs[obs_test] = obs[obs_test]+12
	model_test= model < 6
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
    
def convert_phase_units_actual_all(obs,periods):			
    for i in range(len(obs)):
    
        if obs[i] < 0:
            obs[i] = np.pi+ (np.pi-np.abs(obs[i]))
        elif obs[i] > 0:
            obs[i] = obs[i]

        con_num = periods[i]

        convert_factor = con_num/(2*np.pi)
        obs[i] = obs[i]*convert_factor
	
    return obs
	
def correct_select_daily_phase_actual_single(obs,lon_c,lat_c,lon_e,lat_e, obs_lon):
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
    
def detrend(data,degree):
    detrended=data.copy()
    for i in range(degree,len(data)-degree):
        detrended[i] = data[i] - np.mean(data[i-degree:i+degree])
    return detrended
			
def log_bin_mean(limit, n,x_array,y_array):
    #to get exact numbers areas of len n, must add 1 to n. 
    #eg 10 log spaced numbers returns 9 areas to smooth between. Must add 1 to get 10 areas
    n+=1
    
    #generate log spaced numbers
    result = [1]
    if n>1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result)<n:
        next_value = result[-1]*ratio
        if next_value - result[-1] >= 1:
            # safe zone. next_value will be a different integer
            result.append(next_value)
        else:
            # problem! same integer. we need to find next_value by artificially incrementing previous value
            result.append(result[-1]+1)
            # recalculate the ratio so that the remaining values will scale correctly
            ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    # round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
    log_spaced =  np.array(map(lambda x: round(x)-1, result),dtype=np.uint64)
    #amend last value to make it size of input array
    print log_spaced
    log_spaced = log_spaced.tolist()
    log_spaced[-1] = np.uint64(limit)
    
    #start smoothing of arrays
    smoothed_x_array = np.array([])
    smoothed_y_array = np.array([])
    
    x_array = x_array[::-1]
    y_array = y_array[::-1]
    
    for i in range(len(log_spaced)-1):
        try:
            start_i+=1
            end_i+=1
            print np.mean(x_array[log_spaced[start_i]:log_spaced[end_i]])
            smoothed_x_array = np.append(smoothed_x_array,np.mean(x_array[log_spaced[start_i]:log_spaced[end_i]]))
            smoothed_y_array = np.append(smoothed_y_array,np.mean(y_array[log_spaced[start_i]:log_spaced[end_i]]))
        except:
            start_i = 0
            end_i = 1   
            print np.mean(x_array[log_spaced[start_i]:log_spaced[end_i]])
            smoothed_x_array = np.append(smoothed_x_array,np.mean(x_array[log_spaced[start_i]:log_spaced[end_i]]))
            smoothed_y_array = np.append(smoothed_y_array,np.mean(y_array[log_spaced[start_i]:log_spaced[end_i]]))
            
    return smoothed_x_array,smoothed_y_array			

def log_bin_median(limit, n,x_array,y_array):
    #to get exact numbers areas of len n, must add 1 to n. 
    #eg 10 log spaced numbers returns 9 areas to smooth between. Must add 1 to get 10 areas
    n+=1
    
    #generate log spaced numbers
    result = [1]
    if n>1:  # just a check to avoid ZeroDivisionError
        ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    while len(result)<n:
        next_value = result[-1]*ratio
        if next_value - result[-1] >= 1:
            # safe zone. next_value will be a different integer
            result.append(next_value)
        else:
            # problem! same integer. we need to find next_value by artificially incrementing previous value
            result.append(result[-1]+1)
            # recalculate the ratio so that the remaining values will scale correctly
            ratio = (float(limit)/result[-1]) ** (1.0/(n-len(result)))
    # round, re-adjust to 0 indexing (i.e. minus 1) and return np.uint64 array
    log_spaced =  np.array(map(lambda x: round(x)-1, result),dtype=np.uint64)
    #amend last value to make it size of input array
    print log_spaced
    log_spaced = log_spaced.tolist()
    log_spaced[-1] = np.uint64(limit)
    
    #start smoothing of arrays
    smoothed_x_array = np.array([])
    smoothed_y_array = np.array([])
    
    for i in range(len(log_spaced)-1):
        try:
            start_i+=1
            end_i+=1
            smoothed_x_array = np.append(smoothed_x_array,np.median(x_array[log_spaced[start_i]:log_spaced[end_i]]))
            smoothed_y_array = np.append(smoothed_y_array,np.median(y_array[log_spaced[start_i]:log_spaced[end_i]]))
        except:
            start_i = 0
            end_i = 1    
            smoothed_x_array = np.append(smoothed_x_array,np.median(x_array[log_spaced[start_i]:log_spaced[end_i]]))
            smoothed_y_array = np.append(smoothed_y_array,np.median(y_array[log_spaced[start_i]:log_spaced[end_i]]))
            
    return smoothed_x_array,smoothed_y_array

def sidelobe_peak_remove(fb,fr,fi,closest,crit_percent,periods):
    critical_val = (fb[closest]/100.) * crit_percent # certain % of top peak
    up_peak_val = fb[closest]
    down_peak_val = fb[closest]
    i=1
    while up_peak_val >= critical_val:
        if i != 1:
            try:
                up_inds.append(closest+(i-1))
            except:
                up_inds = [closest+(i-1)]
    
        up_peak_val = fb[closest+i]
        i+=1 
    if i == 2:
        up_inds = np.array([])
    i=1
    while down_peak_val >= critical_val:
        if i != 1:
            try:
                down_inds.append(closest-(i-1))
            except:
                down_inds = [closest-(i-1)]
        down_peak_val = fb[closest-i]
        i+=1

    if i == 2:
        down_inds = []
    all_inds = np.concatenate((down_inds,up_inds))
    all_inds = np.sort(all_inds)
    all_inds = [int(i) for i in all_inds]
    
    #set significant indices for mag, fr, and fi to 0
    fb[closest] = 0
    fr[closest] = 0
    fi[closest] = 0
    fb[all_inds] = 0
    fr[all_inds] = 0
    fi[all_inds] = 0
    
    print 'Altering periods: ', periods[closest], periods[all_inds] 
    
    return fb,fr,fi
    
def phase_offset_correct(key_period,periods,ph):
    closest_period_index = min(range(len(periods)), key=lambda i: abs(periods[i]-key_period))

    print 'ph = ', ph[closest_period_index-5:closest_period_index+5]
    print 'diffs = ', ph[closest_period_index-5]-ph[closest_period_index-4]
    print 'diffs = ', ph[closest_period_index-3]-ph[closest_period_index-2]
    
    closest_period = periods[closest_period_index]
    closest_phase = ph[closest_period_index]
    print 'closest period = ',closest_period
    print 'closest phase = ',closest_phase
    
    if periods[closest_period_index] == key_period:
        return closest_phase    
    
    elif periods[closest_period_index] > key_period:
        highest_period = periods[closest_period_index]
        lowest_period = periods[closest_period_index+1]
        #print 'highest_period = ',highest_period
        #print 'lowest_period = ',lowest_period
        phase_1 = ph[closest_period_index]
        phase_2 = ph[closest_period_index+1]
        period_diff = highest_period-lowest_period
        small_period_diff = closest_period-key_period
        period_ratio = period_diff/small_period_diff
        period_percent = 100./period_ratio
        
    
    else:
        highest_period = periods[closest_period_index-1]
        lowest_period = periods[closest_period_index]
        #print 'highest_period = ',highest_period
        #print 'lowest_period = ',lowest_period
        phase_1 = ph[closest_period_index]
        phase_2 = ph[closest_period_index-1]
        period_diff = highest_period-lowest_period
        #print 'period diff = ', period_diff
        small_period_diff = key_period-closest_period
        #print 'small period diff = ', small_period_diff
        period_ratio = period_diff/small_period_diff
        #print 'period ratio = ', period_ratio
        period_percent = 100./period_ratio
    
    highest_phase = np.max([phase_1,phase_2])
    lowest_phase = np.min([phase_1,phase_2])

    
    print 'highest phase = ',highest_phase
    print 'lowest phase = ', lowest_phase
    print 'period percent = ',period_percent
    phase_diff_test = highest_phase-lowest_phase
    print highest_phase-lowest_phase
    if phase_diff_test > np.pi:
        phase_diff_1 = (np.pi*2)-highest_phase
        phase_diff = phase_diff_1+lowest_phase
        wrap = True
    else:
        phase_diff = highest_phase-lowest_phase
        wrap = False
    print 'Phase diff = ', phase_diff
    phase_percent = phase_diff/100.
    phase_offset = phase_percent*period_percent
    print 'Phase offset = ', phase_offset
    
    print 'wrap = ', wrap
    if wrap == True:
        if closest_phase == lowest_phase:
            corrected_phase = closest_phase-phase_offset
            if corrected_phase < 0:
                corrected_phase = (2*np.pi)-np.abs(corrected_phase)
        else:
            corrected_phase = closest_phase+phase_offset
            if corrected_phase > (2*np.pi):
                gap = corrected_phase-(2*np.pi)
                corrected_phase = 0 + gap
                

    else:
        if closest_phase == lowest_phase:
            corrected_phase = closest_phase+phase_offset
        else:
            corrected_phase = closest_phase-phase_offset

    print 'Corrected phase = ', corrected_phase

    return corrected_phase

def phase_start_point_correct(period,ph,valid_times):
    full_cycle = np.pi*2

    start_time = valid_times[0]
    
    if start_time >= period:
        while start_time >= period:
            start_time = start_time - period
    fraction_off = start_time/period
    offset = full_cycle*fraction_off
    
    ph = ph + offset

    if ph > np.pi:
        ph = -np.pi + (ph - np.pi)
    
    return ph

def phase_start_point_correct_all(periods,phase,valid_times):
    full_cycle = np.pi*2
    
    for i in range(len(periods)):
        #correct phases for offset from start point
        start_time = valid_times[0]
        if start_time >= periods[i]:
            while start_time >= periods[i]:
                start_time = start_time - periods[i]
        fraction_off = start_time/periods[i]
        offset = full_cycle*fraction_off
    
        phase[i] = phase[i] + offset

        if phase[i] > np.pi:
            phase[i] = -np.pi + (phase[i] - np.pi)
    
    return phase


def phase_start_correct(time):
    time = time-time[0]
    return time
    
def try_int(s):
    "Convert to integer if possible."
    try: return int(s)
    except: return s

def natsort_key(s):
    "Used internally to get a tuple by which s is sorted."
    import re
    return map(try_int, re.findall(r'(\d+|\D+)', s))

def natcmp(a, b):
    "Natural string comparison, case sensitive."
    return cmp(natsort_key(a), natsort_key(b))

def natcasecmp(a, b):
    "Natural string comparison, ignores case."
    return natcmp(a.lower(), b.lower())

def natsort(seq, cmp=natcmp):
    "In-place natural string sort."
    seq.sort(cmp)
    
def natsorted(seq, cmp=natcmp):
    "Returns a copy of seq, sorted by natural string sort."
    import copy
    temp = copy.copy(seq)
    natsort(temp, cmp)
    return temp

def sinc(xs):
    pxs = np.pi*xs
    return np.where(np.fabs(pxs)<1e-3, 1.0-pxs*pxs/6.0, np.sin(pxs)/pxs)

def hanning_window(xs, halfwidth):
    win =  0.5 + 0.5*np.cos(np.pi*xs/halfwidth)
    return np.where(np.fabs(xs)<=halfwidth, win, 0.0)

def periodic_interp(fr,fi,zoomfact,periods,key_period,n,amp_corr,window='hanning', alpha=6.0):
            
    closest_period_index = min(range(len(periods)), key=lambda i: abs(periods[i]-key_period))

    peak_fr = fr[closest_period_index-5:closest_period_index+6]
    peak_fi = fi[closest_period_index-5:closest_period_index+6]
    
    peak_periods = periods[closest_period_index-5:closest_period_index+6]
    
    peak_real = peak_fr/n
    peak_real = peak_real*2
    peak_real = peak_real*amp_corr

    peak_imag = peak_fi/n
    peak_imag = peak_imag*2
    peak_imag = peak_imag*amp_corr
    
    reverse_peak_periods = peak_periods[::-1]
    data = peak_real[::-1]
    data2 = peak_imag[::-1]
    
    #interpolate real component
    zoomfact = int(zoomfact)
    if (zoomfact < 1):
        print "zoomfact must be >= 1."
        return 0.0
    elif zoomfact==1:
        return data
    newN = len(data)*zoomfact
    # Space out the data
    comb = np.zeros((zoomfact, len(data)), dtype='d')
    comb[0] += data
    comb = np.reshape(np.transpose(comb), (newN,))
    # Compute the offsets
    xs = np.zeros(newN, dtype='d')
    xs[:newN/2+1] = np.arange(newN/2+1, dtype='d')/zoomfact
    xs[-newN/2:]  = xs[::-1][newN/2-1:-1]
    # Calculate the sinc times window for the kernel
    if window.lower()=="kaiser":
        win = _window_function[window](xs, len(data)/2, alpha)
    else:
        win = hanning_window(xs, len(data)/2)
    kernel = win * sinc(xs)
    newreal = FFT.irfft(FFT.rfft(kernel) * FFT.rfft(comb))
    
    #interpolate imag component
    if (zoomfact < 1):
        print "zoomfact must be >= 1."
        return 0.0
    elif zoomfact==1:
        return data2
    newN = len(data2)*zoomfact
    # Space out the data
    comb = np.zeros((zoomfact, len(data2)), dtype='d')
    comb[0] += data2
    comb = np.reshape(np.transpose(comb), (newN,))
    # Compute the offsets
    xs = np.zeros(newN, dtype='d')
    xs[:newN/2+1] = np.arange(newN/2+1, dtype='d')/zoomfact
    xs[-newN/2:]  = xs[::-1][newN/2-1:-1]
    # Calculate the sinc times window for the kernel
    if window.lower()=="kaiser":
        win = _window_function[window](xs, len(data2)/2, alpha)
    else:
        win = hanning_window(xs, len(data2)/2)
    kernel = win * sinc(xs)
    newimag = FFT.irfft(FFT.rfft(kernel) * FFT.rfft(comb))
    
    #get new periods for interpolation
    new_x = []
    ns = 0
    ne = 1
    over_its = len(reverse_peak_periods)-1
    for i in range(over_its):
        new_x = np.append(new_x,np.linspace(reverse_peak_periods[ns],reverse_peak_periods[ne],zoomfact,endpoint=False))
        ns+=1
        ne+=1

    new_x = np.append(new_x,reverse_peak_periods[-1])

    #get max mag in selected range (and associated phase)
    newreal = newreal[:len(new_x)]
    newimag = newimag[:len(new_x)]
    newmag = []
    newphase = []
    for i in range(len(newreal)):
        newmag.append(np.abs(complex(newreal[i],newimag[i])))
        newphase.append(np.angle(complex(newreal[i],newimag[i])))
    peakind = np.argmax(newmag)

    amp = newmag[peakind]
    phase = newphase[peakind]
    
    return amp,phase
    
def haversine(lat1, lon1, lat2, lon2):
    """
    Calculate the great circle distance between two points 
    on the earth (specified in decimal degrees)
    """
    # convert decimal degrees to radians 
    lon1, lat1, lon2, lat2 = map(radians, [lon1, lat1, lon2, lat2]) 
    # haversine formula 
    dlon = lon2 - lon1  
    dlat = lat2 - lat1  
    a = sin(dlat/2)**2 + cos(lat1) * cos(lat2) * sin(dlon/2)**2
    c = 2 * asin(sqrt(a)) 
    km = 6367 * c
    return km
    
def orthogonal_1_1_res_var(x_array,y_array): # x3,y3 is the point
    
    x = np.arange(0,1000,1)
    y = np.arange(0,1000,1)
    
    x1 = x[0] 
    y1 = y[0]
    x2 = x[-1]
    y2 = y[-1]
    
    all_dists = []
    
    for i in range(len(x_array)):
        x3 = x_array[i] 
        y3 = y_array[i]
    
        px = x2-x1
        py = y2-y1

        something = px*px + py*py

        u =  ((x3 - x1) * px + (y3 - y1) * py) / float(something)

        if u > 1:
            u = 1
        elif u < 0:
            u = 0

        x = x1 + u * px
        y = y1 + u * py

        dx = x - x3
        dy = y - y3

        # Note: If the actual distance does not matter,
        # if you only want to compare what this function
        # returns to other results of this function, you
        # can just return the squared distance instead
        # (i.e. remove the sqrt) to gain a little performance

        dist = math.sqrt(dx*dx + dy*dy)
        all_dists = np.append(all_dists,dist)
    
    dist_sum = np.sum(all_dists)
    dist_ave = np.average(all_dists)
        
    return np.around(dist_sum,2),np.around(dist_ave,2)
    

def anthrome_classify(anthload,obs_lat,obs_lon):
    
    anthromes = {
    11: 'Urban',
    12: 'Mixed settlements',
    21: 'Rice villages',
    22: 'Irrigated villages',
    23: 'Rainfed villages',
    24: 'Pastoral villages',
    31: 'Residential irrigated croplands',
    32: 'Residential rainfed croplands',
    33: 'Populated croplands',
    34: 'Remote croplands',
    41: 'Residential rangelands',
    42: 'Populated rangelands',
    43: 'Remote rangelands',
    51: 'Residential woodlands',
    52: 'Populated woodlands',
    53: 'Remote woodlands',
    54: 'Inhabited treeless and barren lands',
    61: 'Wild woodlands',
    62: 'Wild treeless and barren lands'
    }

    #Anthrome classes may be grouped for analysis into Anthrome levels:
    simple_anthromes = {
    11:	'Dense Settlements',
    12:	'Dense Settlements',
    21:	'Villages',
    22:	'Villages',
    23:	'Villages',
    24:	'Villages',
    31:	'Croplands',
    32:	'Croplands',
    33:	'Croplands',
    34:	'Croplands',
    41:	'Rangelands',
    42:	'Rangelands',
    43:	'Rangelands',
    51:	'Seminatural',
    52:	'Seminatural',
    53:	'Seminatural',
    54:	'Seminatural',
    61:	'Wildlands',
    62:	'Wildlands',
    }

    simple_groups = {
    1: 'Dense Settlement',
    2: 'Village',
    3: 'Cropland',
    4: 'Rangeland',
    5: 'Seminaturnal',
    6: 'Wildland',
    }

    #anthload = Dataset(anthfile)
    
    anth_lon_c = anthload.variables['lon'][:]
    anth_lat_c = anthload.variables['lat'][:]
    anthdata = anthload.variables['anthro2_a2000'][:]
    
    anth_lat_c = anth_lat_c[::-1]
    
    anthdata = anthdata[::-1,:]
    anthdata = anthdata//10
    
    lon_diff = (anth_lon_c[1]-anth_lon_c[0])/2.
    lat_diff = (anth_lat_c[1] - anth_lat_c[0])/2.
    
    anth_lon_e = anth_lon_c-lon_diff
    anth_lon_e = np.append(anth_lon_e,180.00)
    anth_lat_e = anth_lat_c-lat_diff
    anth_lat_e = np.append(anth_lat_e,90.00)

    if len(obs_lat) > 1:
        validity = []
        class_name = []

        for o_lat,o_lon in zip(obs_lat,obs_lon): 
            lat_i,lon_i = obs_model_gridbox(anth_lat_e,anth_lon_e,o_lat,o_lon) 
            anthdata_grid = anthdata[lat_i,lon_i]
            
            try:
                class_n = simple_groups[anthdata_grid]
            except:
                anthdata_grid = 0
                class_n = 'Ocean'
    
            if anthdata_grid != 1:
                valid = 'valid'
            else:
                valid = 'invalid'
                
            validity.append(valid)
            class_name.append(class_n)
        validity = np.array(validity)
        class_name = np.array(class_name)
    
    
    else:
        obs_lon = obs_lon[0]
        obs_lat = obs_lat[0]
        lat_i,lon_i = obs_model_gridbox(anth_lat_e,anth_lon_e,obs_lat,obs_lon)

        anthdata_grid = anthdata[lat_i,lon_i]
        #print obs_lon,obs_lat
        #print anthdata_grid

        #fig, ax = plt.subplots(1,1)
        #m = Basemap(projection='cyl',llcrnrlat=-90,urcrnrlat=90,\
        #                llcrnrlon=-180,\
        #                 urcrnrlon=180,\
        #                resolution='h')
                        
        #m.drawcoastlines(linewidth=0.5)
        #m.drawcountries()
        #x,y=m(anth_lon_c,anth_lat_c)
        #cmap = mpl.cm.get_cmap('rainbow', 7)    # 7 discrete colors
        #norm = mpl.colors.BoundaryNorm(np.arange(8),7)
        #plt.pcolormesh(x,y,anthdata, cmap=cmap, norm=norm)
        #plt.colorbar()
        #plt.show()
    
        try:
            print 'Site is %s'%(simple_groups[anthdata_grid])
            class_name = simple_groups[anthdata_grid]
        except:
            print 'Site is Ocean'
            anthdata_grid = 0
            class_name = 'Ocean'
    
        if anthdata_grid != 1:
            validity = 'valid'
        else:
            validity = 'invalid'
        
    
    return validity,class_name
    
def find_nearest_point_index(x_array,y_array,x_val,y_val):
    point_diff = []
    for i in range(len(x_array)):
        x_diff = np.abs(x_array[i] - x_val)
        y_diff = np.abs(y_array[i] - y_val)
        
        point_diff.append(np.sqrt((x_diff**2)+(y_diff**2)))
        
    return np.argmin(point_diff)   
    
def get_country(lat, lon):
    data = json.load(urllib2.urlopen('http://maps.googleapis.com/maps/api/geocode/json?latlng=%s,%s&sensor=false' % (lat, lon)))
    if len(data['results']) >= 1:
        for result in data['results']:
            for component in result['address_components']:
                if 'country' in component['types']:
                    return component['long_name']        
    else:
        return 'Ocean'

def get_country_2(lat, lon):
    location = ''
    
    while location == '':
        geolocator = Nominatim()
        location = geolocator.reverse((lat,lon))
        location_dict = location.raw
    
    try:
        address_dict = location_dict['address']
        country_code = address_dict['country_code']
        country = transformations.cc_to_cn('ly')
        return country
        
    except:
        return 'Ocean'        

        
def get_continent(name):
    if name == 'Ocean':
        return 'Ocean'
    else:
        try:
            continent = transformations.cn_to_ctn(name)
        except:
            continent = name
        return continent
    
    
def detect_peaks(x, mph=None, mpd=1, threshold=0, edge='rising',
                 kpsh=False, valley=False, show=False, ax=None):

    """Detect peaks in data based on their amplitude and other features.

    Parameters
    ----------
    x : 1D array_like
        data.
    mph : {None, number}, optional (default = None)
        detect peaks that are greater than minimum peak height.
    mpd : positive integer, optional (default = 1)
        detect peaks that are at least separated by minimum peak distance (in
        number of data).
    threshold : positive number, optional (default = 0)
        detect peaks (valleys) that are greater (smaller) than `threshold`
        in relation to their immediate neighbors.
    edge : {None, 'rising', 'falling', 'both'}, optional (default = 'rising')
        for a flat peak, keep only the rising edge ('rising'), only the
        falling edge ('falling'), both edges ('both'), or don't detect a
        flat peak (None).
    kpsh : bool, optional (default = False)
        keep peaks with same height even if they are closer than `mpd`.
    valley : bool, optional (default = False)
        if True (1), detect valleys (local minima) instead of peaks.
    show : bool, optional (default = False)
        if True (1), plot data in matplotlib figure.
    ax : a matplotlib.axes.Axes instance, optional (default = None).

    Returns
    -------
    ind : 1D array_like
        indeces of the peaks in `x`.

    Notes
    -----
    The detection of valleys instead of peaks is performed internally by simply
    negating the data: `ind_valleys = detect_peaks(-x)`
    
    The function can handle NaN's 

    See this IPython Notebook [1]_.

    References
    ----------
    .. [1] http://nbviewer.ipython.org/github/demotu/BMC/blob/master/notebooks/DetectPeaks.ipynb

    Examples
    --------
    >>> from detect_peaks import detect_peaks
    >>> x = np.random.randn(100)
    >>> x[60:81] = np.nan
    >>> # detect all peaks and plot data
    >>> ind = detect_peaks(x, show=True)
    >>> print(ind)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # set minimum peak height = 0 and minimum peak distance = 20
    >>> detect_peaks(x, mph=0, mpd=20, show=True)

    >>> x = [0, 1, 0, 2, 0, 3, 0, 2, 0, 1, 0]
    >>> # set minimum peak distance = 2
    >>> detect_peaks(x, mpd=2, show=True)

    >>> x = np.sin(2*np.pi*5*np.linspace(0, 1, 200)) + np.random.randn(200)/5
    >>> # detection of valleys instead of peaks
    >>> detect_peaks(x, mph=0, mpd=20, valley=True, show=True)

    >>> x = [0, 1, 1, 0, 1, 1, 0]
    >>> # detect both edges
    >>> detect_peaks(x, edge='both', show=True)

    >>> x = [-2, 1, -2, 2, 1, 1, 3, 0]
    >>> # set threshold = 2
    >>> detect_peaks(x, threshold = 2, show=True)
    """

    x = np.atleast_1d(x).astype('float64')
    if x.size < 3:
        return np.array([], dtype=int)
    if valley:
        x = -x
    # find indices of all peaks
    dx = x[1:] - x[:-1]
    # handle NaN's
    indnan = np.where(np.isnan(x))[0]
    if indnan.size:
        x[indnan] = np.inf
        dx[np.where(np.isnan(dx))[0]] = np.inf
    ine, ire, ife = np.array([[], [], []], dtype=int)
    if not edge:
        ine = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) > 0))[0]
    else:
        if edge.lower() in ['rising', 'both']:
            ire = np.where((np.hstack((dx, 0)) <= 0) & (np.hstack((0, dx)) > 0))[0]
        if edge.lower() in ['falling', 'both']:
            ife = np.where((np.hstack((dx, 0)) < 0) & (np.hstack((0, dx)) >= 0))[0]
    ind = np.unique(np.hstack((ine, ire, ife)))
    # handle NaN's
    if ind.size and indnan.size:
        # NaN's and values close to NaN's cannot be peaks
        ind = ind[np.in1d(ind, np.unique(np.hstack((indnan, indnan-1, indnan+1))), invert=True)]
    # first and last values of x cannot be peaks
    if ind.size and ind[0] == 0:
        ind = ind[1:]
    if ind.size and ind[-1] == x.size-1:
        ind = ind[:-1]
    # remove peaks < minimum peak height
    if ind.size and mph is not None:
        ind = ind[x[ind] >= mph]
    # remove peaks - neighbors < threshold
    if ind.size and threshold > 0:
        dx = np.min(np.vstack([x[ind]-x[ind-1], x[ind]-x[ind+1]]), axis=0)
        ind = np.delete(ind, np.where(dx < threshold)[0])
    # detect small peaks closer than minimum peak distance
    if ind.size and mpd > 1:
        ind = ind[np.argsort(x[ind])][::-1]  # sort ind by peak height
        idel = np.zeros(ind.size, dtype=bool)
        for i in range(ind.size):
            if not idel[i]:
                # keep peaks with the same height if kpsh is True
                idel = idel | (ind >= ind[i] - mpd) & (ind <= ind[i] + mpd) \
                    & (x[ind[i]] > x[ind] if kpsh else True)
                idel[i] = 0  # Keep current peak
        # remove the small peaks and sort back the indices by their occurrence
        ind = np.sort(ind[~idel])

    if show:
        if indnan.size:
            x[indnan] = np.nan
        if valley:
            x = -x
        _plot(x, mph, mpd, threshold, edge, valley, ax, ind)

    return ind
    
def moving_average(x, n, type='exponential'):
    x = np.asarray(x)
    if type=='simple':
        weights = np.ones(n)
    else:
        weights = np.exp(np.linspace(-1., 0., n))

    weights /= weights.sum()

    a =  np.convolve(x, weights, mode='full')[:len(x)]
    a[:n] = a[n]
    return a
    
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')

def spectra_breakpoints_fixed(periods,mag,ofac):

    bp1 = 1.
    bp2 = 10.

    #cut periods/mag to key range to get spectral slope estimation
    end_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-0.))
    start_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-100.))
    bp_periods = periods[start_cut_i:end_cut_i+1]
    bp_mag = mag[start_cut_i:end_cut_i+1]

    #remove points n around quarter diurnal peak
    diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-0.25))
    n_points = int((320*ofac)+np.floor(ofac/2.))
    rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
    #bp_periods = np.delete(bp_periods,rm_points)
    #bp_mag = np.delete(bp_mag,rm_points)
    bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    #remove points n around third diurnal peak
    diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-0.33333333))
    n_points = int((240*ofac)+np.floor(ofac/2.))
    rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
    #bp_periods = np.delete(bp_periods,rm_points)
    #bp_mag = np.delete(bp_mag,rm_points)
    bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    #remove points n around half diurnal peak
    diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-0.5))
    n_points = int((160*ofac)+np.floor(ofac/2.))
    rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
    #bp_periods = np.delete(bp_periods,rm_points)
    #bp_mag = np.delete(bp_mag,rm_points)
    bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    #remove points n around diurnal peak
    diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1.))
    n_points = int((100*ofac)+np.floor(ofac/2.))
    rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
    #bp_periods = np.delete(bp_periods,rm_points)
    #bp_mag = np.delete(bp_mag,rm_points)
    bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    #make periods (and mag) go in ascending order.
    bp_periods = bp_periods[::-1]
    bp_mag = bp_mag[::-1]
    
    #get ave power of each regime
    test = bp_periods < bp1
    bp_periods1 = bp_periods[test]
    bp_mag1 = bp_mag[test]
    power1 = np.average(bp_mag1)
    len1 = len(bp_periods1)
    test = (bp_periods >= bp1) & (bp_periods < bp2)
    bp_periods2 = bp_periods[test]
    bp_mag2 = bp_mag[test]
    power2 = np.average(bp_mag2)
    len2 = len(bp_periods2)
    test = bp_periods >= bp2
    bp_periods3 = bp_periods[test]
    len3 = len(bp_periods3)
    bp_mag3 = bp_mag[test]
    power3 = np.average(bp_mag3)
    
    #log periods and mag, will keep log distributed nature after interpolation
    bp_periods = np.log(bp_periods)
    bp_mag = np.log(bp_mag)

    #smooth mag using univariate splines
    s = UnivariateSpline(bp_periods, bp_mag, k = 1, s=(6000*ofac))
    xs = np.linspace(bp_periods[0], bp_periods[-1],300)
    ys = s(xs)
    bp_periods = np.copy(xs)
    bp_mag = np.copy(ys)
    
    #remove first point as may be spurious
    bp_periods = bp_periods[1:]
    bp_mag = bp_mag[1:]

    bp1 = np.log(1.)
    bp2 = np.log(10.)

    test = bp_periods < bp1
    bp_periods1 = bp_periods[test]
    bp_mag1 = bp_mag[test]
    len1 = len(bp_periods1)
    test = (bp_periods >= bp1) & (bp_periods < bp2)
    bp_periods2 = bp_periods[test]
    bp_mag2 = bp_mag[test]
    len2 = len(bp_periods2)
    test = bp_periods >= bp2
    bp_periods3 = bp_periods[test]
    len3 = len(bp_periods3)
    bp_mag3 = bp_mag[test]

    x_lens = [len1,len2,len3]

    popt, pcov = curve_fit(spectra_f,x_lens,bp_mag)
    
    grad1 = popt[0]
    grad2 = popt[1]
    grad3 = popt[2]
    yinter = popt[3]
    
    return grad1,grad2,grad3,power1,power2,power3
    
def getdof(iwin):
    #c50 relates to type of window used, 0 is rectangular, 2 is hanning
    c50 = [0.500, 0.344, 0.167, 0.250, 0.096]

    c2 = 2.0 * c50[iwin] * c50[iwin]
    denom = 1.0 + c2 - c2/ 1.
    neff = 1. / denom
    dof = 2.0 * neff
    return dof
    
    
def getchi2(dof,alpha):
    tol = 1.0e-3
    itmax = 100
    ierr = 0
    
#use approximation for dof > 30 (Eq. 1.132 in Sachs (1984))
#----------------------------------------------------------
    if dof > 30.0:
        #za = -getz(alpha)   # NB: Eq. requires change of sign for percentile
        if ierr == 1:
            return
        x = 2.0 / 9.0 / dof
        chi2 = dof * (1.0 - x + za * np.sqrt(x))**3.0
    else:
        iter = 0
        lm = 0.0
        rm = 1000.0
        if alpha > 0.5:
            eps = (1.0 - alpha) * tol
        else:
            eps = alpha * tol
        while iter <= itmax:
            iter= iter + 1
            chi2 = 0.5 * (lm + rm)
            ac = 1.0 - scipy.special.gammainc(0.5*dof, 0.5*chi2)
            if np.abs(ac - alpha) <= eps:
                break
            if ac > alpha:
                lm = chi2
            else:
                rm = chi2

    if iter > itmax:
        print "Error in GETCHI2: Iter > ItMax"
        ierr = 1
        return
    
    getchi2 = chi2

    return getchi2

def chi_signif(periods,window_t_or_f):
    if window_t_or_f == 't':
        num = 2
    else:
        num = 0
    
    #get degrees of freedom
    dof = getdof(num) # 0 is rectangular,2 is hanning window

    #get scaling factors for red noise model
    fac80 = getchi2(dof, 0.20) / dof
    fac85 = getchi2(dof, 0.15) / dof
    fac90 = getchi2(dof, 0.10) / dof
    fac95 = getchi2(dof, 0.05) / dof
    fac99 = getchi2(dof, 0.01) / dof
    fac99_9 = getchi2(dof, 0.001) / dof
    fac99_99 = getchi2(dof, 0.0001) / dof
    
    # critical false alarm level after Thomson (1990)
    # -----------------------------------------------
    ave_seg_len = len(periods)
    alphacrit = 1.0 / ave_seg_len
    faccrit = getchi2(dof, alphacrit) / dof
    
    return fac80,fac85,fac90,fac95,fac99,fac99_9,fac99_99,faccrit
    
    
def spectra_f_3(x,grad1,grad2,grad3,y_inter):

    x1_r = x[0]
    x2_r = x[1]
    x3_r = x[2]  
   
    diff1 = x2_r[0] - x1_r[-1]
    diff2 = x3_r[0] - x2_r[-1]
   
    x1_r = x1_r - x1_r[0]
    x2_r = x2_r - x2_r[0]
    x3_r = x3_r - x3_r[0]
   
    p1 = grad1*x1_r+y_inter
    y_inter2 = grad1*(x1_r[-1]+diff1)+y_inter
    p2 = grad2*x2_r+y_inter2
    y_inter3 = grad2*(x2_r[-1]+diff2)+y_inter2
    p3 = grad3*x3_r+y_inter3

    p = [p1]+[p2]+[p3]
    p = [item for sublist in p for item in sublist]
    
    return p
    
def spectra_f_2(x,grad1,grad2,y_inter):

    x1_r = x[0]
    x2_r = x[1]
   
    diff1 = x2_r[0] - x1_r[-1]
   
    x1_r = x1_r - x1_r[0]
    x2_r = x2_r - x2_r[0]
   
    p1 = grad1*x1_r+y_inter
    y_inter2 = grad1*(x1_r[-1]+diff1)+y_inter
    p2 = grad2*x2_r+y_inter2

    p = [p1]+[p2]
    p = [item for sublist in p for item in sublist]
    
    return p    

def spectra_f_parts_3(x,grad1,grad2,grad3,y_inter):
    
    x1_r = x[0]
    x2_r = x[1]
    x3_r = x[2]  
   
    diff1 = x2_r[0] - x1_r[-1]
    diff2 = x3_r[0] - x2_r[-1]
   
    x1_r = x1_r - x1_r[0]
    x2_r = x2_r - x2_r[0]
    x3_r = x3_r - x3_r[0]
   
    p1 = grad1*x1_r+y_inter
    y_inter2 = grad1*(x1_r[-1]+diff1)+y_inter
    p2 = grad2*x2_r+y_inter2
    y_inter3 = grad2*(x2_r[-1]+diff2)+y_inter2
    p3 = grad3*x3_r+y_inter3
  
    return p1,p2,p3
    
def spectra_f_parts_2(x,grad1,grad2,y_inter):
    
    x1_r = x[0]
    x2_r = x[1]
   
    diff1 = x2_r[0] - x1_r[-1]
   
    x1_r = x1_r - x1_r[0]
    x2_r = x2_r - x2_r[0]
   
    p1 = grad1*x1_r+y_inter
    y_inter2 = grad1*(x1_r[-1]+diff1)+y_inter
    p2 = grad2*x2_r+y_inter2
  
    return p1,p2

    
def spectra_fit_fixed_piecewise(periods,mag,ofac,lower_period_limit,upper_period_limit,breakpoint):

    start_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-upper_period_limit))
    end_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-lower_period_limit))
    
    if periods[end_cut_i] > lower_period_limit:
        end_cut_i = end_cut_i + 1.
    bp_periods = periods[start_cut_i:end_cut_i]
    bp_mag = mag[start_cut_i:end_cut_i]
    
    if (lower_period_limit < 1./10.) & (1./10. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./10.))
        n_points = int((500*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./9.) & (1./9. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./9.))
        n_points = int((450*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./8.) & (1./8. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./8.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./7.) & (1./7. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./7.))
        n_points = int((350*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./6.) & (1./6. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./6.))
        n_points = int((300*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./5.) & (1./5. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./5.))
        n_points = int((250*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./4.) & (1./4. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./4.))
        n_points = int((200*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./3.) & (1./3. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./3.))
        n_points = int((150*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./2.) & (1./2. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./2.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit < 1.) & (1. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 365.25/4.) & (365.25/4. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/4.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit < 365.25/3.) & (365.25/3. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/3.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 365.25/2.) & (365.25/2. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/2.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
     
    if (lower_period_limit < 365.25) & (365.25 < upper_period_limit):
        annual_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(annual_i-n_points,(annual_i+n_points)+1)
        if rm_points[0] < 0:
            rm_points = range(0,(annual_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    #make periods (and mag) go in ascending order.
    bp_periods = bp_periods[::-1]
    bp_mag = bp_mag[::-1]
    
    #log periods and mag, will keep log distributed nature after interpolation
    bp_periods = np.log(bp_periods)
    bp_mag = np.log(bp_mag)
    
    bp1 = np.log(breakpoint)
    
    test = bp_periods < bp1

    p1 = bp_periods[test]
    bp_mag1 = bp_mag[test]
    len1 = len(p1)
    test = bp_periods >= bp1
    p2 = bp_periods[test]
    bp_mag2 = bp_mag[test]
    len2 = len(p2)
    
    x_lens = [p1,p2]
    
    popt, pcov = curve_fit(spectra_f_2,x_lens,bp_mag)
    
    grad1 = popt[0]
    grad2 = popt[1]
    yinter = popt[2]

    full_periods = np.log(periods[::-1])
    test = full_periods < bp1
    p1 = full_periods[test]
    test = full_periods >= bp1
    p2 = full_periods[test]
    x_lens = [p1,p2]
    
    calc_1,calc_2 = spectra_f_parts_2(x_lens,grad1,grad2,yinter)
    
    calc_1 = calc_1.tolist()
    calc_2 = calc_2.tolist()
    
    calc_y = calc_1+calc_2 
    
    #calculate stats of fitted lines
    #-----------------------------------------
  
    ave1 = np.average(calc_1)
    ave2 = np.average(calc_2)
    med1 = np.median(calc_1)
    med2 = np.median(calc_2)
    sum1 = np.sum(calc_1)
    sum2 = np.sum(calc_2)
    
    y1_s = calc_1[0]
    y1_e = calc_1[-1]
    y2_s = calc_2[0]
    y2_e = calc_2[-1]
    
    plt.loglog(periods,mag)
    plt.loglog(np.exp(p1),np.exp(calc_1))
    plt.loglog(np.exp(p2),np.exp(calc_2))
    plt.show()
    
    return grad1,grad2,np.exp(bp1),np.exp(p1),np.exp(calc_1),np.exp(p2),np.exp(calc_2),np.exp(ave1),np.exp(ave2),np.exp(med1),np.exp(med2),np.exp(sum1),np.exp(sum2),np.exp(y1_s),np.exp(y1_e),np.exp(y2_s),np.exp(y2_e)

def spectra_fit_fixed_linear(periods,mag,ofac,lower_period_limit,upper_period_limit,breakpoint):

    start_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-upper_period_limit))
    end_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-lower_period_limit))
    
    if periods[end_cut_i] > lower_period_limit:
        end_cut_i = end_cut_i + 1.
    bp_periods = periods[start_cut_i:end_cut_i]
    bp_mag = mag[start_cut_i:end_cut_i]
    
    if (lower_period_limit < 1./10.) & (1./10. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./10.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./9.) & (1./9. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./9.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./8.) & (1./8. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./8.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./7.) & (1./7. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./7.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./6.) & (1./6. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./6.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./5.) & (1./5. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./5.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./4.) & (1./4. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./4.))
        n_points = int((320*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./3.) & (1./3. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./3.))
        n_points = int((240*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./2.) & (1./2. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./2.))
        n_points = int((160*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit >= 1.) & (1. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit >= 365.25/4.) & (365.25/4. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/4.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit >= 365.25/3.) & (365.25/3. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/3.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit >= 365.25/2.) & (365.25/2. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/2.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
     
    if (lower_period_limit >= 365.25) & (365.25 < upper_period_limit):
        annual_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(annual_i-n_points,(annual_i+n_points)+1)
        if rm_points[0] < 0:
            rm_points = range(0,(annual_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    #make periods (and mag) go in ascending order.
    bp_periods = bp_periods[::-1]
    bp_mag = bp_mag[::-1]
    
    #log periods and mag, will keep log distributed nature after interpolation
    bp_periods = np.log(bp_periods)
    bp_mag = np.log(bp_mag)
    
    bp1 = np.log(breakpoint)

    test = bp_periods < bp1
    p1 = bp_periods[test]
    bp_mag1 = bp_mag[test]
    test = bp_periods >= bp1
    p2 = bp_periods[test]
    bp_mag2 = bp_mag[test]
    
    grad1, yinter1, r_value, p_value, std_err = stats.linregress(p1,bp_mag1)
    grad2, yinter2, r_value, p_value, std_err = stats.linregress(p2,bp_mag2)
    
    full_periods = np.log(periods[::-1])
    test = full_periods < bp1
    p1 = full_periods[test]
    test = full_periods >= bp1
    p2 = full_periods[test]   
    
    calc_1 = grad1*p1+yinter1
    calc_2 = grad2*p2+yinter2
    
    calc_1 = calc_1.tolist()
    calc_2 = calc_2.tolist()
    
    calc_y = calc_1+calc_2 
    
    orig_periods_log = np.log(periods)[::-1]
    
    #calculate stats of fitted lines
    #-----------------------------------------
  
    ave1 = np.average(calc_1)
    ave2 = np.average(calc_2)
    med1 = np.median(calc_1)
    med2 = np.median(calc_2)
    sum1 = np.sum(calc_1)
    sum2 = np.sum(calc_2)
    
    y1_s = calc_1[0]
    y1_e = calc_1[-1]
    y2_s = calc_2[0]
    y2_e = calc_2[-1]
    
    return grad1,grad2,np.exp(bp1),np.exp(p1),np.exp(calc_1),np.exp(p2),np.exp(calc_2),np.exp(ave1),np.exp(ave2),np.exp(med1),np.exp(med2),np.exp(sum1),np.exp(sum2),np.exp(y1_s),np.exp(y1_e),np.exp(y2_s),np.exp(y2_e)

def spectra_fit_free_piecewise(periods,mag,ofac,lower_period_limit,upper_period_limit,lower_break,upper_break):
    
    start_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-upper_period_limit))
    end_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-lower_period_limit))
    
    if periods[end_cut_i] > lower_period_limit:
        end_cut_i = end_cut_i + 1.
    bp_periods = periods[start_cut_i:end_cut_i]
    bp_mag = mag[start_cut_i:end_cut_i]
    
    if (lower_period_limit < 1./10.) & (1./10. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./10.))
        n_points = int((500*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./9.) & (1./9. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./9.))
        n_points = int((450*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./8.) & (1./8. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./8.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./7.) & (1./7. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./7.))
        n_points = int((350*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./6.) & (1./6. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./6.))
        n_points = int((300*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./5.) & (1./5. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./5.))
        n_points = int((250*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./4.) & (1./4. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./4.))
        n_points = int((200*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./3.) & (1./3. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./3.))
        n_points = int((150*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./2.) & (1./2. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./2.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit < 1.) & (1. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 365.25/4.) & (365.25/4. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/4.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit < 365.25/3.) & (365.25/3. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/3.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 365.25/2.) & (365.25/2. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/2.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
     
    if (lower_period_limit < 365.25) & (365.25 < upper_period_limit):
        annual_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(annual_i-n_points,(annual_i+n_points)+1)
        if rm_points[0] < 0:
            rm_points = range(0,(annual_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    #make periods (and mag) go in ascending order.
    bp_periods = bp_periods[::-1]
    bp_mag = bp_mag[::-1]
    
    #log periods and mag, will keep log distributed nature after interpolation
    bp_periods = np.log(bp_periods)
    bp_mag = np.log(bp_mag)
    
    #set limiting ranges for breakpoints
    breakpoint_test = (bp_periods >= np.log(lower_break)) & (bp_periods < np.log(upper_break))
    bp_periods_1 = bp_periods[breakpoint_test]
    bp1 = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-bp_periods_1[0]))

    all_error1 = []
    all_error2 = []
    all_bp1 = []
    all_grad1 = []
    all_grad2 = []
    all_yinter = []

    for i in range(len(bp_periods_1)):
        bp1 = bp_periods_1[i]
    
        test = bp_periods < bp1
        bp_periods1 = bp_periods[test]
        bp_mag1 = bp_mag[test]
        len1 = len(bp_periods1)
        test = bp_periods >= bp1
        bp_periods2 = bp_periods[test]
        bp_mag2 = bp_mag[test]
        len2 = len(bp_periods2)
    
        x_lens = [bp_periods1,bp_periods2]
    
        popt, pcov = curve_fit(spectra_f_2,x_lens,bp_mag)
    
        grad1 = popt[0]
        grad2 = popt[1]
        yinter = popt[2]
    
        bp_mag_est = spectra_f_2(x_lens,grad1,grad2,yinter)    
        
        line1,line2 = spectra_f_parts_2(x_lens,grad1,grad2,yinter)
        
        err1=np.sqrt(np.sum((np.abs(line1-bp_mag1))**2)/len(bp_mag1))
        err2=np.sqrt(np.sum((np.abs(line2-bp_mag2))**2)/len(bp_mag2))
    
        all_error1.append(err1)
        all_error2.append(err2)
        all_bp1.append(np.average((bp_periods1[-1],bp_periods2[0])))
        all_grad1.append(grad1)
        all_grad2.append(grad2)
        all_yinter.append(yinter)
        
    
    #get best fit based on best rank of errors for each line
    all_error1 = np.array(all_error1)
    all_error2 = np.array(all_error2)

    temp = all_error1.argsort()
    error1_ranks = np.empty(len(all_error1), int)
    error1_ranks[temp] = np.arange(len(all_error1))
    temp = all_error2.argsort()
    error2_ranks = np.empty(len(all_error2), int)
    error2_ranks[temp] = np.arange(len(all_error2))

    bp_joint_error_rank = error1_ranks+error2_ranks
    bp_ind_err = np.argmin(bp_joint_error_rank)

    bp1 = all_bp1[bp_ind_err]
    grad1 = all_grad1[bp_ind_err]
    grad2 = all_grad2[bp_ind_err]
    yinter = all_yinter[bp_ind_err]
    
    full_periods = np.log(periods[::-1])
    test = full_periods < bp1
    p1 = full_periods[test]
    test = full_periods >= bp1
    p2 = full_periods[test]    
    
    x_lens = [p1,p2]
    
    calc_1,calc_2 = spectra_f_parts_2(x_lens,grad1,grad2,yinter)

    calc_1 = calc_1.tolist()
    calc_2 = calc_2.tolist()
    
    calc_y = calc_1+calc_2 
    
    orig_periods_log = np.log(periods)[::-1]
    
    #calculate stats of fitted lines
    #-----------------------------------------
    ave1 = np.average(calc_1)
    ave2 = np.average(calc_2)
    med1 = np.median(calc_1)
    med2 = np.median(calc_2)
    sum1 = np.sum(calc_1)
    sum2 = np.sum(calc_2)
    
    y1_s = calc_1[0]
    y1_e = calc_1[-1]
    y2_s = calc_2[0]
    y2_e = calc_2[-1]

    return grad1,grad2,np.exp(bp1),np.exp(p1),np.exp(calc_1),np.exp(p2),np.exp(calc_2),np.exp(ave1),np.exp(ave2),np.exp(med1),np.exp(med2),np.exp(sum1),np.exp(sum2),np.exp(y1_s),np.exp(y1_e),np.exp(y2_s),np.exp(y2_e)

def spectra_fit_free_linear(periods,mag,ofac,lower_period_limit,upper_period_limit,lower_break,upper_break):
    
    start_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-upper_period_limit))
    end_cut_i = min(range(len(periods)), key=lambda i: abs(periods[i]-lower_period_limit))
    
    if periods[end_cut_i] > lower_period_limit:
        end_cut_i = end_cut_i + 1.
    bp_periods = periods[start_cut_i:end_cut_i]
    bp_mag = mag[start_cut_i:end_cut_i]
    
    if (lower_period_limit < 1./10.) & (1./10. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./10.))
        n_points = int((500*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./9.) & (1./9. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./9.))
        n_points = int((450*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./8.) & (1./8. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./8.))
        n_points = int((400*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./7.) & (1./7. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./7.))
        n_points = int((350*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./6.) & (1./6. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./6.))
        n_points = int((300*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./5.) & (1./5. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./5.))
        n_points = int((250*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./4.) & (1./4. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./4.))
        n_points = int((200*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./3.) & (1./3. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./3.))
        n_points = int((150*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 1./2.) & (1./2. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1./2.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit < 1.) & (1. < upper_period_limit):
        diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1.))
        n_points = int((100*ofac)+np.floor(ofac/2.))
        rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 365.25/4.) & (365.25/4. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/4.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    if (lower_period_limit < 365.25/3.) & (365.25/3. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/3.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    if (lower_period_limit < 365.25/2.) & (365.25/2. < upper_period_limit):
        ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25/2.))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
     
    if (lower_period_limit < 365.25) & (365.25 < upper_period_limit):
        annual_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25))
        n_points = int((1*ofac)+np.floor(ofac/2.))
        rm_points = range(annual_i-n_points,(annual_i+n_points)+1)
        if rm_points[0] < 0:
            rm_points = range(0,(annual_i+n_points)+1)
        bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    #make periods (and mag) go in ascending order.
    bp_periods = bp_periods[::-1]
    bp_mag = bp_mag[::-1]
    
    #log periods and mag, will keep log distributed nature after interpolation
    bp_periods = np.log(bp_periods)
    bp_mag = np.log(bp_mag)
    
    #set limiting ranges for breakpoints
    breakpoint_test = (bp_periods >= np.log(lower_break)) & (bp_periods < np.log(upper_break))
    bp_periods_1 = bp_periods[breakpoint_test]
    bp1 = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-bp_periods_1[0]))

    all_error1 = []
    all_error2 = []
    all_bp1 = []
    all_grad1 = []
    all_grad2 = []
    all_yinter = []
    all_yinter2 = []

    for i in range(len(bp_periods_1)):
        bp1 = bp_periods_1[i]
    
        test = bp_periods < bp1
        bp_periods1 = bp_periods[test]
        bp_mag1 = bp_mag[test]
        len1 = len(bp_periods1)
        test = bp_periods >= bp1
        bp_periods2 = bp_periods[test]
        bp_mag2 = bp_mag[test]
        len2 = len(bp_periods2)
    
        grad1, yinter1, r_value, p_value, std_err = stats.linregress(p1,bp_mag1)
        grad2, yinter2, r_value, p_value, std_err = stats.linregress(p2,bp_mag2)
    
        line1 = grad1*bp_periods1+yinter1
        line2 = grad2*bp_periods2+yinter2
        
        err1=np.sqrt(np.sum((np.abs(line1-bp_mag1))**2)/len(bp_mag1))
        err2=np.sqrt(np.sum((np.abs(line2-bp_mag2))**2)/len(bp_mag2))
    
        all_error1.append(err1)
        all_error2.append(err2)
        all_bp1.append(np.average((bp_periods1[-1],bp_periods2[0])))
        all_grad1.append(grad1)
        all_grad2.append(grad2)
        all_yinter.append(yinter)
        all_yinter2.append(yinter2)
        
    #get best fit based on best rank of errors for each line
    all_error1 = np.array(all_error1)
    all_error2 = np.array(all_error2)

    temp = all_error1.argsort()
    error1_ranks = np.empty(len(all_error1), int)
    error1_ranks[temp] = np.arange(len(all_error1))
    temp = all_error2.argsort()
    error2_ranks = np.empty(len(all_error2), int)
    error2_ranks[temp] = np.arange(len(all_error2))

    bp_joint_error_rank = error1_ranks+error2_ranks
    bp_ind_err = np.argmin(bp_joint_error_rank)

    bp1 = all_bp1[bp_ind_err]
    grad1 = all_grad1[bp_ind_err]
    grad2 = all_grad2[bp_ind_err]
    yinter = all_yinter[bp_ind_err]
    yinter2 = all_yinter2[bp_ind_err]
    
    full_periods = np.log(periods[::-1])
    test = full_periods < bp1
    p1 = full_periods[test]
    test = full_periods >= bp1
    p2 = full_periods[test]    
    
    calc_1 = grad1*p1+yinter1
    calc_2 = grad2*p2+yinter2
    
    calc_1 = calc_1.tolist()
    calc_2 = calc_2.tolist()
    
    calc_y = calc_1+calc_2 
    
    orig_periods_log = np.log(periods)[::-1]
    
    #calculate stats of fitted lines
    #-----------------------------------------
    ave1 = np.average(calc_1)
    ave2 = np.average(calc_2)
    med1 = np.median(calc_1)
    med2 = np.median(calc_2)
    sum1 = np.sum(calc_1)
    sum2 = np.sum(calc_2)
    
    y1_s = calc_1[0]
    y1_e = calc_1[-1]
    y2_s = calc_2[0]
    y2_e = calc_2[-1]

    return grad1,grad2,np.exp(bp1),np.exp(p1),np.exp(calc_1),np.exp(p2),np.exp(calc_2),np.exp(ave1),np.exp(ave2),np.exp(med1),np.exp(med2),np.exp(sum1),np.exp(sum2),np.exp(y1_s),np.exp(y1_e),np.exp(y2_s),np.exp(y2_e)

def period_convolution(key_periods,time,mag,ph,mean):
    pi2 = np.pi*2
    
    if np.max(key_periods) == 1.:
        max_time = 1.
    elif np.max(key_periods) == 365.25:
        max_time = 365.25
    
    for i in range(len(key_periods)):
        try:
            waveform = waveform + (mag[i]*(np.cos((pi2*time/key_periods[i])-(ph[i]))))
        except:
            waveform = mag[i]*(np.cos((pi2*time/key_periods[i])-(ph[i])))
    
    ff = np.sqrt(np.mean(np.square(waveform)))/np.mean(np.abs(waveform))
    waveform = mean+waveform
    waveform_min = np.min(waveform)
    waveform_max = np.max(waveform)
    waveform_min_ind = np.argmin(waveform)
    waveform_max_ind = np.argmax(waveform)
    mag = (waveform_max-waveform_min)/2.
    ph_min = time[waveform_min_ind]
    ph_max = time[waveform_max_ind]
    ratio = (pi2)/max_time
    ph_min = ratio*ph_min
    ph_max = ratio*ph_max
    
    return mag,ph_min,ph_max,waveform,ff 

def solar_time_phase_corrector(current_phase,hours,time_diff):
    pi2 = np.pi*2.
    ph_ratio = pi2/hours
    ph_offset = time_diff * ph_ratio
    
    new_phase = current_phase + ph_offset
    remainder = np.mod(np.abs(new_phase),pi2)
    
    if new_phase >= pi2:
        phase = remainder
    elif new_phase < 0:
        phase = pi2 - remainder
    else:
        phase = new_phase
    
    return phase
        
def phase_to_radians(period,ph):
    pi2 = np.pi*2.
    ratio = pi2/period
    ph = ph*ratio
   
    return ph
    
def find_closest(A, target):
    #A must be sorted
    idx = A.searchsorted(target)
    idx = np.clip(idx, 1, len(A)-1)
    left = A[idx-1]
    right = A[idx]
    idx -= target - left < right - target
    return idx
    
def remove_periodic(bp_periods,bp_mag,ofac):
    diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-0.25))
    n_points = int((320*ofac)+np.floor(ofac/2.))
    rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
    #bp_periods = np.delete(bp_periods,rm_points)
    #bp_mag = np.delete(bp_mag,rm_points)
    bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    #remove points n around third diurnal peak
    diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-0.33333333))
    n_points = int((240*ofac)+np.floor(ofac/2.))
    rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
    #bp_periods = np.delete(bp_periods,rm_points)
    #bp_mag = np.delete(bp_mag,rm_points)
    bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    #remove points n around half diurnal peak
    diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-0.5))
    n_points = int((160*ofac)+np.floor(ofac/2.))
    rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
    #bp_periods = np.delete(bp_periods,rm_points)
    #bp_mag = np.delete(bp_mag,rm_points)
    bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))

    #remove points n around diurnal peak
    diurnal_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-1.))
    n_points = int((100*ofac)+np.floor(ofac/2.))
    rm_points = range(diurnal_i-n_points,(diurnal_i+n_points)+1)
    #bp_periods = np.delete(bp_periods,rm_points)
    #bp_mag = np.delete(bp_mag,rm_points)
    bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    #remove points n around half-annual peak
    ha_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-182.625))
    n_points = int((1*ofac)+np.floor(ofac/2.))
    rm_points = range(ha_i-n_points,(ha_i+n_points)+1)
    #bp_periods = np.delete(bp_periods,rm_points)
    #bp_mag = np.delete(bp_mag,rm_points)
    bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
     
    #remove points n around annual peak
    annual_i = min(range(len(bp_periods)), key=lambda i: abs(bp_periods[i]-365.25))
    n_points = int((1*ofac)+np.floor(ofac/2.))
    rm_points = range(annual_i-n_points,(annual_i+n_points)+1)
    if rm_points[0] < 0:
        rm_points = range(0,(annual_i+n_points)+1)
    #bp_periods = np.delete(bp_periods,rm_points)
    #bp_mag = np.delete(bp_mag,rm_points)
    bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    
    return bp_periods,bp_mag
    
def regrid_interp(var,orig_lat_e,orig_lat_c,orig_lon_e,orig_lon_c,min_lat,min_lon,data_type,max_ph):

    new_lat_e = np.arange(-90,91,min_lat)
    new_lat_c = np.arange((-90+(min_lat/2.)),90,min_lat)
    new_lon_e = np.arange(-180,181,min_lon)
    new_lon_c = np.arange((-180+(min_lon/2.)),180,min_lon)
    n_boxes = len(new_lat_c)*len(new_lon_c)

    z = np.empty((len(new_lat_c),len(new_lon_c)))

    #find edge of boxes for new lat and new lon edges to left of gridbox new lat and lon edges are in. 
    for lat_i in range(len(new_lat_c)):
        for lon_i in range(len(new_lon_c)): 
        
            exception_lat = False
            exception_lon = False
            
            #print 'new lon edges', new_lon_e[lon_i],new_lon_e[lon_i+1]  
            #print 'new lat edges', new_lat_e[lat_i],new_lat_e[lat_i+1] 
        
            if new_lat_e[lat_i] > orig_lat_e[-1]:
                diff = new_lat_e[lat_i] - orig_lat_e[-1]
                corrected_new_lat_e = orig_lat_e[0] + diff
                lat_box_left_i = np.where(orig_lat_e == np.max(orig_lat_e[orig_lat_e <= corrected_new_lat_e]))[0][0]
            else:
                lat_box_left_i = np.where(orig_lat_e == np.max(orig_lat_e[orig_lat_e <= new_lat_e[lat_i]]))[0][0]
            
            if new_lat_e[lat_i+1] > orig_lat_e[-1]:
                diff = new_lat_e[lat_i+1] - orig_lat_e[-1]
                corrected_new_lat_e = orig_lat_e[0] + diff
                lat_box_right_i = np.where(orig_lat_e == np.min(orig_lat_e[orig_lat_e >= corrected_new_lat_e]))[0][0]
            else:
                lat_box_right_i = np.where(orig_lat_e == np.min(orig_lat_e[orig_lat_e >= new_lat_e[lat_i+1]]))[0][0]
            if lat_box_right_i - lat_box_left_i < 0: 
                #print 'lat exception'
                exception_lat = True
            
            if new_lon_e[lon_i] > orig_lon_e[-1]:
                diff = new_lon_e[lon_i] - orig_lon_e[-1]
                corrected_new_lon_e = orig_lon_e[0] + diff
                lon_box_left_i = np.where(orig_lon_e == np.max(orig_lon_e[orig_lon_e <= corrected_new_lon_e]))[0][0]   
            else:
                lon_box_left_i = np.where(orig_lon_e == np.max(orig_lon_e[orig_lon_e <= new_lon_e[lon_i]]))[0][0]
            
            if new_lon_e[lon_i+1] > orig_lon_e[-1]:   
                diff = new_lon_e[lon_i+1] - orig_lon_e[-1]
                corrected_new_lon_e = orig_lon_e[0] + diff
                lon_box_right_i = np.where(orig_lon_e == np.min(orig_lon_e[orig_lon_e >= corrected_new_lon_e]))[0][0]
            else:
                lon_box_right_i = np.where(orig_lon_e == np.min(orig_lon_e[orig_lon_e >= new_lon_e[lon_i+1]]))[0][0]
            if lon_box_right_i - lon_box_left_i < 0: 
                #print 'lon exception'
                exception_lon = True
            
            #print 'orig lon edges',orig_lon_e[lon_box_left_i],orig_lon_e[lon_box_right_i]
            #print 'orig lat edges',orig_lat_e[lat_box_left_i],orig_lat_e[lat_box_right_i]
            
            orig_lat_left = orig_lat_e[lat_box_left_i]
            orig_lon_left = orig_lon_e[lon_box_left_i]
            
            new_lat_left = new_lat_e[lat_i]
            new_lat_right = new_lat_e[lat_i+1]
            new_lon_left = new_lon_e[lon_i]
            new_lon_right = new_lon_e[lon_i+1]
            
            if exception_lat == True:
                lat_diff_i = len(orig_lat_e) - lat_box_left_i
            else:
                orig_lat_right = orig_lat_e[lat_box_right_i]
                lat_diff_i = lat_box_right_i - lat_box_left_i
            
            if exception_lon == True:    
                lon_diff_i = len(orig_lon_e) - lon_box_left_i
            else: 
                orig_lon_right = orig_lon_e[lon_box_right_i]
                lon_diff_i = lon_box_right_i - lon_box_left_i
                
            #print lat_box_left_i, lat_box_right_i,lon_box_left_i, lon_box_right_i 
            #print lat_diff_i,lon_diff_i
            
            #how many boxes does new boxes overlap?
            #1 box
            if (lat_diff_i == 1) & (lon_diff_i == 1):
                #print '1 box'
            
                if exception_lat == True:
                    orig_box1_lat_diff = 180. -  np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_right_i]])[0])
                else:
                    orig_box1_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_right_i]])[0])
                
                lat_box1_n = lat_box_left_i

                if exception_lon == True:
                    orig_box1_lon_diff = 360. - np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_right_i]])[0])
                else:
                    orig_box1_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_right_i]])[0])
                
                if lon_box_left_i == len(orig_lon_c):
                    lon_box1_n = 0
                else:
                    lon_box1_n = lon_box_left_i
        
                z[lat_i,lon_i] = var[lat_box1_n,lon_box1_n]
                
                #print var[lat_box1_n,lon_box1_n]
                #print z[lat_i,lon_i]
    
            #2 boxes
            elif (lat_diff_i) + (lon_diff_i) == 3:
                #print '2 boxes'
                #2 lat boxes, 1 lon box
                if (lat_diff_i) == 2:
                    if exception_lat == True:
                        orig_box1_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_left_i+1]])[0])
                        orig_box2_lat_diff = 180. - np.abs(np.diff([orig_lat_e[lat_box_left_i+1] , orig_lat_e[lat_box_right_i]])[0])
                    else:
                        orig_box1_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_left_i+1]])[0])
                        orig_box2_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i+1] , orig_lat_e[lat_box_right_i]])[0])
                    
                    lat_box1_n = lat_box_left_i
                    lat_box2_n = lat_box_right_i-1
                
                    if exception_lon == True:
                        orig_box1_lon_diff = 360. - np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_right_i]])[0])
                    else:
                        orig_box1_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_right_i]])[0])
                        

                    if lon_box_left_i == len(orig_lon_c):
                        lon_box1_n = 0
                    else:
                        lon_box1_n = lon_box_left_i
                
                    actual_box1_lat_diff = np.abs(np.diff([new_lat_e[lat_i],orig_lat_e[lat_box_left_i+1]])[0])
                    actual_box2_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i+1],new_lat_e[lat_i+1]])[0])
                    actual_box1_lon_diff = np.abs(np.diff([new_lon_e[lon_i],new_lon_e[lon_i+1]])[0])
                
                    box1_area = orig_box1_lat_diff * orig_box1_lon_diff
                    actual_box1_area = actual_box1_lat_diff * actual_box1_lon_diff
                    box1_fract = (1./box1_area)*actual_box1_area
                    #print 'actual box 1 fract = ',box1_fract
            
                    box2_area = orig_box2_lat_diff * orig_box1_lon_diff
                    actual_box2_area = actual_box2_lat_diff * actual_box1_lon_diff
                    box2_fract = (1./box2_area)*actual_box2_area
                    #print 'actual box 2 fract = ',box2_fract
                    
                    if data_type == 'ph':
                        data = [var[lat_box1_n,lon_box1_n],var[lat_box2_n,lon_box1_n]]
                        weights=[box1_fract,box2_fract]
                        z[lat_i,lon_i] = cylic_weighted_average(data,weights,max_ph)
                    else:
                        z[lat_i,lon_i] = np.average([var[lat_box1_n,lon_box1_n],
                                                     var[lat_box2_n,lon_box1_n]],weights=[box1_fract,box2_fract])
                       
                    #print var[lat_box1_n,lon_box1_n],var[lat_box2_n,lon_box1_n]                          
                    #print z[lat_i,lon_i]
        
                #2 lon boxes, 1 lat box
                elif (lon_diff_i) == 2:
                    if exception_lat == True:
                        orig_box1_lat_diff = 180. - np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_right_i]])[0])
                    else:
                        orig_box1_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_right_i]])[0])
                    
                    lat_box1_n = lat_box_left_i
                
                    if exception_lon == True:
                        orig_box1_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_left_i+1]])[0])
                        orig_box2_lon_diff = 360. - np.abs(np.diff([orig_lon_e[lon_box_left_i+1] , orig_lon_e[lon_box_right_i]])[0])
                    else:
                        orig_box1_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_left_i+1]])[0])
                        orig_box2_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i+1] , orig_lon_e[lon_box_right_i]])[0])
                    
                    if lon_box_left_i == len(orig_lon_c):
                        lon_box1_n = 0
                        lon_box2_n = 1
                    elif lon_box_left_i+1 == len(orig_lon_c):
                        lon_box1_n = lon_box_left_i
                        lon_box2_n = 0
                    else:
                        lon_box1_n = lon_box_left_i
                        lon_box2_n = lon_box_left_i+1 
                    
                    actual_box1_lat_diff = np.abs(np.diff([new_lat_e[lat_i],new_lat_e[lat_i+1]])[0])
                    actual_box1_lon_diff = np.abs(np.diff([new_lon_e[lon_i],orig_lon_e[lon_box_left_i+1]])[0])
                    actual_box2_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i+1],new_lon_e[lon_i+1]])[0])
                
                    box1_area = orig_box1_lat_diff * orig_box1_lon_diff
                    actual_box1_area = actual_box1_lat_diff * actual_box1_lon_diff
                    box1_fract = (1./box1_area)*actual_box1_area
                    #print 'actual box 1 fract = ',box1_fract
            
                    box2_area = orig_box1_lat_diff * orig_box2_lon_diff
                    actual_box2_area = actual_box1_lat_diff * actual_box2_lon_diff
                    box2_fract = (1./box2_area)*actual_box2_area 
                    #print 'actual box 2 fract = ',box2_fract
            
                    if data_type == 'ph':
                        data = [var[lat_box1_n,lon_box1_n],var[lat_box1_n,lon_box2_n]]
                        weights=[box1_fract,box2_fract]
                        z[lat_i,lon_i] = cylic_weighted_average(data,weights,max_ph)
                    else:
                        z[lat_i,lon_i] = np.average([var[lat_box1_n,lon_box1_n],
                                                     var[lat_box1_n,lon_box2_n]],weights=[box1_fract,box2_fract])
                    
                    #print var[lat_box1_n,lon_box1_n],var[lat_box1_n,lon_box2_n]
                    #print z[lat_i,lon_i]
                    
    
            #4 boxes
            elif (lat_diff_i) + (lon_diff_i) == 4:    
            
                #print '4 boxes'
                if exception_lat == True:
                    orig_box1_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_left_i+1]])[0])
                    orig_box2_lat_diff = 180. - np.abs(np.diff([orig_lat_e[lat_box_left_i+1] , orig_lat_e[lat_box_right_i]])[0])
                else:
                    orig_box1_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i] , orig_lat_e[lat_box_left_i+1]])[0])
                    orig_box2_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i+1] , orig_lat_e[lat_box_right_i]])[0])
                
                lat_box1_n = lat_box_left_i
                lat_box2_n = lat_box_right_i-1
            
                if exception_lon == True:
                    orig_box1_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_left_i+1]])[0])
                    orig_box2_lon_diff = 360. - np.abs(np.diff([orig_lon_e[lon_box_left_i+1] , orig_lon_e[lon_box_right_i]])[0])
                else:
                    orig_box1_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i] , orig_lon_e[lon_box_left_i+1]])[0])
                    orig_box2_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i+1] , orig_lon_e[lon_box_right_i]])[0])
                #print orig_box2_lon_diff
                
                if lon_box_left_i == len(orig_lon_c):
                    lon_box1_n = 0
                    lon_box2_n = 1
                elif lon_box_left_i+1 == len(orig_lon_c):
                    lon_box1_n = lon_box_left_i
                    lon_box2_n = 0
                else:
                    lon_box1_n = lon_box_left_i
                    lon_box2_n = lon_box_left_i+1
                
                actual_box1_lat_diff = np.abs(np.diff([new_lat_e[lat_i],orig_lat_e[lat_box_left_i+1]])[0])
                actual_box2_lat_diff = np.abs(np.diff([orig_lat_e[lat_box_left_i+1],new_lat_e[lat_i+1]])[0])
                actual_box1_lon_diff = np.abs(np.diff([new_lon_e[lon_i],orig_lon_e[lon_box_left_i+1]])[0])
                actual_box2_lon_diff = np.abs(np.diff([orig_lon_e[lon_box_left_i+1],new_lon_e[lon_i+1]])[0])
                
                #lower left box of grid
                box1_area = orig_box1_lat_diff * orig_box1_lon_diff
                actual_box1_area = actual_box1_lat_diff * actual_box1_lon_diff
                box1_fract = (1./box1_area)*actual_box1_area
                #print 'actual box 1 fract = ',box1_fract
    
                #lower right box of grid
                box2_area = orig_box1_lat_diff * orig_box2_lon_diff
                actual_box2_area = actual_box1_lat_diff * actual_box2_lon_diff
                box2_fract = (1./box2_area)*actual_box2_area
                #print 'actual box 2 fract = ',box2_fract
                
                #upper left box of grid
                box3_area = orig_box2_lat_diff * orig_box1_lon_diff
                actual_box3_area = actual_box2_lat_diff * actual_box1_lon_diff
                box3_fract = (1./box3_area)*actual_box3_area
                #print 'actual box 3 fract = ',box3_fract
        
                #upper right box of grid
                box4_area = orig_box2_lat_diff * orig_box2_lon_diff
                actual_box4_area = actual_box2_lat_diff * actual_box2_lon_diff
                box4_fract = (1./box4_area)*actual_box4_area
                #print 'actual box 4 fract = ',box4_fract
                
                if data_type == 'ph':
                    data = [var[lat_box1_n,lon_box1_n],var[lat_box1_n,lon_box2_n],var[lat_box2_n,lon_box1_n],var[lat_box2_n,lon_box2_n]]
                    weights=[box1_fract,box2_fract,box3_fract,box4_fract]
                    z[lat_i,lon_i] = cylic_weighted_average(data,weights,max_ph)
                else:
                    z[lat_i,lon_i] = np.average([var[lat_box1_n,lon_box1_n],
                                                 var[lat_box1_n,lon_box2_n],
                                                 var[lat_box2_n,lon_box1_n],
                                                 var[lat_box2_n,lon_box2_n]],weights=[box1_fract,box2_fract,box3_fract,box4_fract])
                                  
                #print var[lat_box1_n,lon_box1_n],var[lat_box1_n,lon_box2_n],var[lat_box2_n,lon_box1_n],var[lat_box2_n,lon_box2_n]
                #print z[lat_i,lon_i]
        
            #print ''                                      
    
    
    return new_lat_e,new_lat_c,new_lon_e,new_lon_c,z
    
def get_obs_info(present_dir):

    paths = present_dir.split("/")
    species = paths[-4]
    years = paths[-2]
    start_year = years[:4]
    end_year = years[5:9]
    model_path = paths[-1] 
    
    data_split = model_path.split('_')
    if len(data_split) == 7:
        model = data_split[0]
        vres = data_split[1]
        version = data_split[2]
        hres = data_split[3]
        met = data_split[4]
        timeres = data_split[5]
        additional = data_split[6]
        fname = '/work/home/db876/observations/surface/%s/process/GLOBAL_%s_%s_%s_%s_%s_PERIODIC.nc'%(species,vres,species,start_year,end_year,timeres)
    
    if len(data_split) == 3:
        vres = data_split[1]
        timeres = data_split[2]
        fname = '/work/home/db876/observations/surface/%s/process/GLOBAL_%s_%s_%s_%s_%s_PERIODIC.nc'%(species,vres,species,start_year,end_year,timeres)
        
    return fname,species,start_year,end_year,vres,timeres
    
def get_model_info(present_dir):
    
    paths = present_dir.split("/")
    species = paths[-4]
    file_years = paths[-3]
    file_start_year = file_years[:4]
    file_end_year = file_years[5:9]
    years = paths[-2]
    start_year = years[:4]
    end_year = years[5:9]
    model_path = paths[-1] 

    data_split = model_path.split('_') 
    model = data_split[0]
    vres = data_split[1]
    version = data_split[2]
    hres = data_split[3]
    met = data_split[4]
    timeres = data_split[5]
    additional = data_split[6]

    fname = '/work/home/db876/plotting_tools/model_files/%s_%s_%s_%s_%s_%s_%s_%s_%s.nc'%(model,vres,file_start_year,file_end_year,version,hres,met,timeres,additional)
    
    if species == 'GMAOTEMP':
        species = 'GMAO_TEMP'
    if species == 'GMAOUWND':
        species = 'GMAO_UWND'
    if species == 'GMAOVWND':
        species = 'GMAO_VWND'
    if species == 'WINDDIRECTION':
        species = 'WIND_DIRECTION'
    if species == 'WINDSPEED':
        species = 'WIND_SPEED'
    
    return fname,species,start_year,end_year

def clicker_interactive_map_obsmodel(event,species,lat_e,lon_e,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_d_waveform,obs_s_waveform,obs_all_waveform,model_d_waveform,model_s_waveform,model_all_waveform,fig,all_m):
    
    #ind of point clicked near
    ind = event.ind
    print event
    print ind
    print obs_refs[ind]
    print len(obs_refs)
    print obs_refs[:20]
    ind = ind[0]
    ref = obs_refs[ind]
    
    
    #read in obs and model periodic parts
    obs_site_group = obs_period_grp.groups[ref]
    obs_d_mag = obs_site_group.daily_amplitude
    obs_d_ph = obs_site_group.daily_phase
    obs_s_mag = obs_site_group.seasonal_amplitude
    obs_s_ph = obs_site_group.seasonal_phase
    obs_ave = obs_site_group.average
    obs_d_waveform = obs_d_waveform[ind]
    obs_s_waveform = obs_s_waveform[ind]
    obs_all_waveform = obs_all_waveform[ind]
    
    model_site_group = model_period_grp.groups[ref]
    model_d_mag = model_site_group.daily_amplitude
    model_d_ph = model_site_group.daily_phase
    model_s_mag = model_site_group.seasonal_amplitude
    model_s_ph = model_site_group.seasonal_phase
    model_ave = model_site_group.average
    model_d_waveform = model_d_waveform[ind]
    model_s_waveform = model_s_waveform[ind]
    model_all_waveform = model_all_waveform[ind]
      
    #get obs timeseries for site clicked and obs details
    obs_site_group = obs_ts_grp.groups[ref]
    obs_var = obs_site_group.variables[species.lower()][:]
    group = obs_site_group.process_group
    lat = obs_site_group.latitude
    lon = obs_site_group.longitude
    lon = obs_site_group.longitude
    alt = obs_site_group.altitude
    complete = obs_site_group.data_completeness
    a_class = obs_site_group.anthrome_class
    r_class = obs_site_group.raw_class
    continent = loc_dict[tags[ind]]
    country = obs_site_group.country
    obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
    obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
    obs_var_pd = pd.Series(obs_var_mask, index=obs_time_pd)
    
    print lat,lon
    
    #highlight points on click
    pl = all_m[0].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20)
    pl2 = all_m[1].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20)
    pl3 = all_m[2].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20) 
    pl4 = all_m[3].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20) 
    
    #get model timeseries for site clicked
    model_var = model_ts_grp.variables[species.lower()][:]*1e9
    lat_n,lon_n = obs_model_gridbox(lat_e,lon_e,lat,lon)
    model_var_pick = model_var[:,lat_n,lon_n]
    model_var_mask = np.ma.masked_where(model_var_pick<=0,model_var_pick)
    model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
    model_var_pd = pd.Series(model_var_mask, index=model_time_pd)     
    
    #get cut times and datetimes
    d_datetime = obs_datetimes[:24]
    s_datetime = obs_datetimes[:8766]
    d_time_pd = pd.date_range(start = d_datetime[0],end = d_datetime[-1], freq = 'H')
    s_time_pd = pd.date_range(start = s_datetime[0],end = s_datetime[-1], freq = 'H')
    
    #make panda series
    obs_d_waveform_pd = pd.Series(obs_d_waveform, index=d_time_pd)
    obs_s_waveform_pd = pd.Series(obs_s_waveform, index=s_time_pd)
    obs_all_waveform_pd = pd.Series(obs_all_waveform, index=obs_time_pd)
    
    model_d_waveform_pd = pd.Series(model_d_waveform, index=d_time_pd)
    model_s_waveform_pd = pd.Series(model_s_waveform, index=s_time_pd)
    model_all_waveform_pd = pd.Series(model_all_waveform, index=model_time_pd)
    
    fig.canvas.draw()
        
    fig2 = plt.figure(figsize =(24,12))
    fig2.patch.set_facecolor('white')
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=1)
    ax3 = plt.subplot2grid((2,2), (1,1), rowspan=1)

    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_var_pd, color='black',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Raw Data')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='red',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Raw Data')
    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_all_waveform_pd, color='gray',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Significant Waveform')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_all_waveform_pd, color='lightcoral',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Significant Waveform')
    ax2.plot_date(d_time_pd.to_pydatetime(), obs_d_waveform_pd, color='black',marker = 'o', markersize = 6,markeredgecolor = None, mew = 0, linestyle = '--', label = 'Obs. Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), obs_s_waveform_pd, color='black',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Obs. Seasonal Waveform') 
    ax2.plot_date(d_time_pd.to_pydatetime(), model_d_waveform_pd, color='red',marker = 'o', markersize = 6, markeredgecolor = None, mew = 0,linestyle = '--', label = 'Model Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), model_s_waveform_pd, color='red',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Model Seasonal Waveform')   
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=model_d_waveform_pd>obs_d_waveform_pd, facecolor='yellow', interpolate=True)
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=obs_d_waveform_pd>model_d_waveform_pd, facecolor='blue', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=model_s_waveform_pd>obs_s_waveform_pd, facecolor='yellow', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=obs_s_waveform_pd>model_s_waveform_pd, facecolor='blue', interpolate=True)
    
    ax1.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s, Grid Index = %s,%s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,lat_n,lon_n))
    
    #plt.legend(loc = 'lower right') 
    plt.tight_layout()
    ax1.grid()
    ax2.grid()
    ax3.grid()
    
    pl.pop(0).remove()
    pl2.pop(0).remove()
    pl3.pop(0).remove()
    pl4.pop(0).remove()
    
    plt.show()
    
    return

def clicker_interactive_map_model(event,plot_type,lat_e,lon_e,linear_lats,linear_lons,datetimes,all_periods,var,d_waveform,s_waveform,all_waveform,fig,m):
    
    #ind of point clicked near
    ind = event.ind
    ind = ind[0]
    
    #highlight points on click
    pl = m.plot([linear_lons[ind]], [linear_lats[ind]], 's', ms=20, alpha=0.6, color='yellow',zorder=20)
    if plot_type == 'multi':
        pl2 = 1
        pl3 = 1
        pl4 = 1
    
    #get model timeseries for site clicked
    lat_n,lon_n = obs_model_gridbox(lat_e,lon_e,linear_lats[ind],linear_lons[ind])
    var_pick = var[:,lat_n,lon_n]   
    var_mask = np.ma.masked_where(var_pick<=0,var_pick)
    time_pd = pd.date_range(start = datetimes[0],end = datetimes[-1], freq = 'H')
    var_pd = pd.Series(var_mask, index=time_pd) 
    
    #read in periodic parts
    d_mag = all_periods[4,lat_n,lon_n]
    d_ph = all_periods[9,lat_n,lon_n]
    s_mag = all_periods[14,lat_n,lon_n]
    s_ph = all_periods[19,lat_n,lon_n]
    ave = all_periods[20,lat_n,lon_n]
    d_waveform = d_waveform[ind]
    s_waveform = s_waveform[ind]
    all_waveform = all_waveform[ind]
    
    #get cut times and datetimes
    d_datetime = datetimes[:24]
    s_datetime = datetimes[:8766]
    d_time_pd = pd.date_range(start = d_datetime[0],end = d_datetime[-1], freq = 'H')
    s_time_pd = pd.date_range(start = s_datetime[0],end = s_datetime[-1], freq = 'H')
    
    #make panda series
    d_waveform_pd = pd.Series(d_waveform, index=d_time_pd)
    s_waveform_pd = pd.Series(s_waveform, index=s_time_pd)
    all_waveform_pd = pd.Series(all_waveform, index=time_pd)
    
    fig.canvas.draw()
        
    fig2 = plt.figure(figsize =(24,12))
    fig2.patch.set_facecolor('white')
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=1)
    ax3 = plt.subplot2grid((2,2), (1,1), rowspan=1)

    ax1.plot_date(time_pd.to_pydatetime(), var_pd, color='red', markersize = 3,marker = 'o', markeredgecolor = None, mew = 0, label = 'Model Raw Data')
    ax1.plot_date(time_pd.to_pydatetime(), all_waveform_pd, color='lightcoral', markersize = 3, marker = 'o', markeredgecolor = None, mew = 0, label = 'Model Significant Waveform')
    ax2.plot_date(d_time_pd.to_pydatetime(), d_waveform_pd, color='red', markersize = 6, marker = 'o', markeredgecolor = None, mew = 0, label = 'Model Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), s_waveform_pd, color='red', markersize = 3, marker = 'o', markeredgecolor = None, mew = 0, label = 'Model Seasonal Waveform') 
    
    #plt.legend(loc = 'lower right') 
    plt.tight_layout()
    ax1.grid()
    ax2.grid()
    ax3.grid()
    
    pl.pop(0).remove()
    
    plt.show()
    
    return
    
def clicker_interactive_xy_obsmodel_single(event,species,lat_e,lon_e,obs_lats,obs_lons,date,time,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_chosen,model_chosen,obs_d_waveform,obs_s_waveform,obs_all_waveform,model_d_waveform,model_s_waveform,model_all_waveform,fig,ax):
    
    #ind of point clicked near
    ind = event.ind
    ind = ind[0]
    ref = obs_refs[ind]
    
    #read in obs and model periodic parts
    obs_site_group = obs_period_grp.groups[ref]
    obs_d_mag = obs_site_group.daily_amplitude
    obs_d_ph = obs_site_group.daily_phase
    obs_s_mag = obs_site_group.seasonal_amplitude
    obs_s_ph = obs_site_group.seasonal_phase
    obs_ave = obs_site_group.average
    obs_d_waveform = obs_d_waveform[ind]
    obs_s_waveform = obs_s_waveform[ind]
    obs_all_waveform = obs_all_waveform[ind]
    obs_d_ff = obs_site_group.daily_ff
    obs_s_ff = obs_site_group.seasonal_ff
    
    model_site_group = model_period_grp.groups[ref]
    model_d_mag = model_site_group.daily_amplitude
    model_d_ph = model_site_group.daily_phase
    model_s_mag = model_site_group.seasonal_amplitude
    model_s_ph = model_site_group.seasonal_phase
    model_ave = model_site_group.average
    model_d_waveform = model_d_waveform[ind]
    model_s_waveform = model_s_waveform[ind]
    model_all_waveform = model_all_waveform[ind]
    model_d_ff = model_site_group.daily_ff
    model_s_ff = model_site_group.seasonal_ff
    
    #highlight points on click
    pl = ax.plot(obs_chosen[ind], model_chosen[ind], 'o', ms=20, alpha=0.6, color='yellow',zorder=20)        
    
    #get obs timeseries for site clicked and obs details
    obs_site_group = obs_ts_grp.groups[ref]
    obs_var = obs_site_group.variables[species.lower()][:]
    group = obs_site_group.process_group
    lat = obs_site_group.latitude
    lon = obs_site_group.longitude
    lon = obs_site_group.longitude
    alt = obs_site_group.altitude
    complete = obs_site_group.data_completeness
    a_class = obs_site_group.anthrome_class
    r_class = obs_site_group.raw_class
    continent = loc_dict[tags[ind]]
    country = obs_site_group.country
    obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
    obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
    obs_var_pd = pd.Series(obs_var_mask, index=obs_time_pd)
    
    #get model timeseries for site clicked
    model_var = model_ts_grp.variables[species.lower()][:]*1e9
    lat_n,lon_n = obs_model_gridbox(lat_e,lon_e,obs_lats[ind],obs_lons[ind])
    model_var_pick = model_var[:,lat_n,lon_n]
    model_var_mask = np.ma.masked_where(model_var_pick<=0,model_var_pick)
    model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
    model_var_pd = pd.Series(model_var_mask, index=model_time_pd)     
    
    #get cut times and datetimes
    d_datetime = obs_datetimes[:24]
    s_datetime = obs_datetimes[:8766]
    d_time_pd = pd.date_range(start = d_datetime[0],end = d_datetime[-1], freq = 'H')
    s_time_pd = pd.date_range(start = s_datetime[0],end = s_datetime[-1], freq = 'H')
    
    #make panda series
    obs_d_waveform_pd = pd.Series(obs_d_waveform, index=d_time_pd)
    obs_s_waveform_pd = pd.Series(obs_s_waveform, index=s_time_pd)
    obs_all_waveform_pd = pd.Series(obs_all_waveform, index=obs_time_pd)
    
    model_d_waveform_pd = pd.Series(model_d_waveform, index=d_time_pd)
    model_s_waveform_pd = pd.Series(model_s_waveform, index=s_time_pd)
    model_all_waveform_pd = pd.Series(model_all_waveform, index=model_time_pd)
    
    fig.canvas.draw()
        
    fig2 = plt.figure(figsize =(24,12))
    fig2.patch.set_facecolor('white')
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=1)
    ax3 = plt.subplot2grid((2,2), (1,1), rowspan=1)

    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_var_pd, color='black',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Raw Data')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='red',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Raw Data')
    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_all_waveform_pd, color='gray',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Significant Waveform')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_all_waveform_pd, color='lightcoral',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Significant Waveform')
    ax2.plot_date(d_time_pd.to_pydatetime(), obs_d_waveform_pd, color='black',marker = 'o', markersize = 6,markeredgecolor = None, mew = 0, linestyle = '--', label = 'Obs. Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), obs_s_waveform_pd, color='black',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Obs. Seasonal Waveform') 
    ax2.plot_date(d_time_pd.to_pydatetime(), model_d_waveform_pd, color='red',marker = 'o', markersize = 6, markeredgecolor = None, mew = 0,linestyle = '--', label = 'Model Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), model_s_waveform_pd, color='red',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Model Seasonal Waveform')   
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=model_d_waveform_pd>obs_d_waveform_pd, facecolor='yellow', interpolate=True)
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=obs_d_waveform_pd>model_d_waveform_pd, facecolor='blue', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=model_s_waveform_pd>obs_s_waveform_pd, facecolor='yellow', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=obs_s_waveform_pd>model_s_waveform_pd, facecolor='blue', interpolate=True)
    
    ax1.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s, Grid Index = %s,%s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,lat_n,lon_n))
    
    #plt.legend(loc = 'lower right') 
    plt.tight_layout()
    ax1.grid()
    ax2.grid()
    ax3.grid()
    
    pl.pop(0).remove()
    
    plt.show()
    
    return

def clicker_interactive_xy_obsmodel_multi(event,species,lat_e,lon_e,obs_lats,obs_lons,date,time,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_d_waveform,obs_s_waveform,obs_all_waveform,model_d_waveform,model_s_waveform,model_all_waveform,fig,ax,ax2,ax3,ax4,ax5,ax7,ax8):
    
    #ind of point clicked near
    ind = event.ind
    ind = ind[0]
    ref = obs_refs[ind]
    
    #read in obs and model periodic parts
    obs_site_group = obs_period_grp.groups[ref]
    obs_d_mag = obs_site_group.daily_amplitude
    obs_d_ph = obs_site_group.daily_phase
    obs_s_mag = obs_site_group.seasonal_amplitude
    obs_s_ph = obs_site_group.seasonal_phase
    obs_ave = obs_site_group.average
    obs_d_ff = obs_site_group.daily_ff
    obs_s_ff = obs_site_group.seasonal_ff
    obs_d_waveform = obs_d_waveform[ind]
    obs_s_waveform = obs_s_waveform[ind]
    obs_all_waveform = obs_all_waveform[ind]
    
    model_site_group = model_period_grp.groups[ref]
    model_d_mag = model_site_group.daily_amplitude
    model_d_ph = model_site_group.daily_phase
    model_s_mag = model_site_group.seasonal_amplitude
    model_s_ph = model_site_group.seasonal_phase
    model_ave = model_site_group.average
    model_d_ff = model_site_group.daily_ff
    model_s_ff = model_site_group.seasonal_ff
    model_d_waveform = model_d_waveform[ind]
    model_s_waveform = model_s_waveform[ind]
    model_all_waveform = model_all_waveform[ind]
    
    
    #highlight points on click
    pl = ax.plot(obs_d_mag, model_d_mag, 'o', ms=12, alpha=0.9, color='black',zorder=20)        
    pl2 = ax2.plot(obs_s_mag, model_s_mag, 'o', ms=12, alpha=0.9, color='black',zorder=20)
    pl3 = ax3.plot(obs_ave, model_ave, 'o', ms=12, alpha=0.9, color='black',zorder=20)
    pl4 = ax4.plot(obs_d_ph, model_d_ph, 'o', ms=12, alpha=0.9, color='black',zorder=20)
    pl5 = ax5.plot(obs_s_ph, model_s_ph, 'o', ms=12, alpha=0.9, color='black',zorder=20)
    pl6 = ax7.plot(obs_d_ff, model_d_ff, 'o', ms=12, alpha=0.9, color='black',zorder=20)
    pl7 = ax8.plot(obs_s_ff, model_s_ff, 'o', ms=12, alpha=0.9, color='black',zorder=20)
    
    #get obs timeseries for site clicked and obs details
    obs_site_group = obs_ts_grp.groups[ref]
    obs_var = obs_site_group.variables[species.lower()][:]
    group = obs_site_group.process_group
    lat = obs_site_group.latitude
    lon = obs_site_group.longitude
    lon = obs_site_group.longitude
    alt = obs_site_group.altitude
    complete = obs_site_group.data_completeness
    a_class = obs_site_group.anthrome_class
    r_class = obs_site_group.raw_class
    continent = loc_dict[tags[ind]]
    country = obs_site_group.country
    obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
    obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
    obs_var_pd = pd.Series(obs_var_mask, index=obs_time_pd)
    
    #get model timeseries for site clicked
    model_var = model_ts_grp.variables[species.lower()][:]*1e9
    lat_n,lon_n = obs_model_gridbox(lat_e,lon_e,obs_lats[ind],obs_lons[ind])
    model_var_pick = model_var[:,lat_n,lon_n]
    model_var_mask = np.ma.masked_where(model_var_pick<=0,model_var_pick)
    model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
    model_var_pd = pd.Series(model_var_mask, index=model_time_pd) 
    
    #get cut times and datetimes
    d_datetime = obs_datetimes[:24]
    s_datetime = obs_datetimes[:8766]
    d_time_pd = pd.date_range(start = d_datetime[0],end = d_datetime[-1], freq = 'H')
    s_time_pd = pd.date_range(start = s_datetime[0],end = s_datetime[-1], freq = 'H')
    
    #make panda series
    obs_d_waveform_pd = pd.Series(obs_d_waveform, index=d_time_pd)
    obs_s_waveform_pd = pd.Series(obs_s_waveform, index=s_time_pd)
    obs_all_waveform_pd = pd.Series(obs_all_waveform, index=obs_time_pd)
    
    model_d_waveform_pd = pd.Series(model_d_waveform, index=d_time_pd)
    model_s_waveform_pd = pd.Series(model_s_waveform, index=s_time_pd)
    model_all_waveform_pd = pd.Series(model_all_waveform, index=model_time_pd)
    
    fig.canvas.draw()
        
    fig2 = plt.figure(figsize =(24,12))
    fig2.patch.set_facecolor('white')
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=1)
    ax3 = plt.subplot2grid((2,2), (1,1), rowspan=1)

    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_var_pd, color='black',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Raw Data')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='red',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Raw Data')
    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_all_waveform_pd, color='gray',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Significant Waveform')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_all_waveform_pd, color='lightcoral',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Significant Waveform')
    ax2.plot_date(d_time_pd.to_pydatetime(), obs_d_waveform_pd, color='black',marker = 'o', markersize = 6,markeredgecolor = None, mew = 0, label = 'Obs. Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), obs_s_waveform_pd, color='black',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Obs. Seasonal Waveform') 
    ax2.plot_date(d_time_pd.to_pydatetime(), model_d_waveform_pd, color='red',marker = 'o', markersize = 6, markeredgecolor = None, mew = 0, label = 'Model Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), model_s_waveform_pd, color='red',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Model Seasonal Waveform')   
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=model_d_waveform_pd>obs_d_waveform_pd, facecolor='yellow', interpolate=True)
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=obs_d_waveform_pd>model_d_waveform_pd, facecolor='blue', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=model_s_waveform_pd>obs_s_waveform_pd, facecolor='yellow', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=obs_s_waveform_pd>model_s_waveform_pd, facecolor='blue', interpolate=True)
    
    ax1.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s, Grid Index = %s,%s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,lat_n,lon_n))
    
    #plt.legend(loc = 'lower right') 
    plt.tight_layout()
    ax1.grid()
    ax2.grid()
    ax3.grid()
    
    pl.pop(0).remove()
    pl2.pop(0).remove()
    pl3.pop(0).remove()
    pl4.pop(0).remove()
    pl5.pop(0).remove()
    
    plt.show()
    
    return
    
def clicker_interactive_map_obsmodel_loc(event,species,lat_e,lon_e,obs_datetimes,model_datetimes,obs_ts_grp,model_ts_grp,obs_period_grp,model_period_grp,obs_refs,tags,loc_dict,obs_d_waveform,obs_s_waveform,obs_all_waveform,model_d_waveform,model_s_waveform,model_all_waveform,fig,all_m):
    
    #ind of point clicked near
    ind = event.ind
    print event
    print ind
    print obs_refs[ind]
    print len(obs_refs)
    print obs_refs[:20]
    ind = ind[0]
    ref = obs_refs[ind]
    
    
    #read in obs and model periodic parts
    obs_site_group = obs_period_grp.groups[ref]
    obs_d_mag = obs_site_group.daily_amplitude
    obs_d_ph = obs_site_group.daily_phase
    obs_s_mag = obs_site_group.seasonal_amplitude
    obs_s_ph = obs_site_group.seasonal_phase
    obs_ave = obs_site_group.average
    obs_d_waveform = obs_d_waveform[ind]
    obs_s_waveform = obs_s_waveform[ind]
    obs_all_waveform = obs_all_waveform[ind]
    
    model_site_group = model_period_grp.groups[ref]
    model_d_mag = model_site_group.daily_amplitude
    model_d_ph = model_site_group.daily_phase
    model_s_mag = model_site_group.seasonal_amplitude
    model_s_ph = model_site_group.seasonal_phase
    model_ave = model_site_group.average
    model_d_waveform = model_d_waveform[ind]
    model_s_waveform = model_s_waveform[ind]
    model_all_waveform = model_all_waveform[ind]
      
    #get obs timeseries for site clicked and obs details
    obs_site_group = obs_ts_grp.groups[ref]
    obs_var = obs_site_group.variables[species.lower()][:]
    group = obs_site_group.process_group
    lat = obs_site_group.latitude
    lon = obs_site_group.longitude
    lon = obs_site_group.longitude
    alt = obs_site_group.altitude
    complete = obs_site_group.data_completeness
    a_class = obs_site_group.anthrome_class
    r_class = obs_site_group.raw_class
    continent = loc_dict[tags[ind]]
    country = obs_site_group.country
    obs_var_mask = np.ma.masked_where(obs_var<=0,obs_var)
    obs_time_pd = pd.date_range(start = obs_datetimes[0],end = obs_datetimes[-1], freq = 'H')
    obs_var_pd = pd.Series(obs_var_mask, index=obs_time_pd)
    
    print lat,lon
    
    #highlight points on click
    pl = all_m[0].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20)
    pl2 = all_m[1].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20)
    pl3 = all_m[2].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20) 
    pl4 = all_m[3].plot(lon, lat, 'o', ms=20, alpha=0.6, color='yellow',zorder=20) 
    
    #get model timeseries for site clicked
    model_var = model_ts_grp.variables[species.lower()][:]*1e9
    lat_n,lon_n = obs_model_gridbox(lat_e,lon_e,lat,lon)
    model_var_pick = model_var[:,lat_n,lon_n]
    model_var_mask = np.ma.masked_where(model_var_pick<=0,model_var_pick)
    model_time_pd = pd.date_range(start = model_datetimes[0],end = model_datetimes[-1], freq = 'H')
    model_var_pd = pd.Series(model_var_mask, index=model_time_pd)     
    
    #get cut times and datetimes
    d_datetime = obs_datetimes[:24]
    s_datetime = obs_datetimes[:8766]
    d_time_pd = pd.date_range(start = d_datetime[0],end = d_datetime[-1], freq = 'H')
    s_time_pd = pd.date_range(start = s_datetime[0],end = s_datetime[-1], freq = 'H')
    
    #make panda series
    obs_d_waveform_pd = pd.Series(obs_d_waveform, index=d_time_pd)
    obs_s_waveform_pd = pd.Series(obs_s_waveform, index=s_time_pd)
    obs_all_waveform_pd = pd.Series(obs_all_waveform, index=obs_time_pd)
    
    model_d_waveform_pd = pd.Series(model_d_waveform, index=d_time_pd)
    model_s_waveform_pd = pd.Series(model_s_waveform, index=s_time_pd)
    model_all_waveform_pd = pd.Series(model_all_waveform, index=model_time_pd)
    
    fig.canvas.draw()
        
    fig2 = plt.figure(figsize =(24,12))
    fig2.patch.set_facecolor('white')
    ax1 = plt.subplot2grid((2,2), (0,0), colspan=2)
    ax2 = plt.subplot2grid((2,2), (1,0), colspan=1)
    ax3 = plt.subplot2grid((2,2), (1,1), rowspan=1)

    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_var_pd, color='black',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Raw Data')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_var_pd, color='red',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Raw Data')
    ax1.plot_date(obs_time_pd.to_pydatetime(), obs_all_waveform_pd, color='gray',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Obs. Significant Waveform')
    ax1.plot_date(model_time_pd.to_pydatetime(), model_all_waveform_pd, color='lightcoral',marker = 'o', markeredgecolor = None, mew = 0, markersize = 3, label = 'Model Significant Waveform')
    ax2.plot_date(d_time_pd.to_pydatetime(), obs_d_waveform_pd, color='black',marker = 'o', markersize = 6,markeredgecolor = None, mew = 0, linestyle = '--', label = 'Obs. Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), obs_s_waveform_pd, color='black',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Obs. Seasonal Waveform') 
    ax2.plot_date(d_time_pd.to_pydatetime(), model_d_waveform_pd, color='red',marker = 'o', markersize = 6, markeredgecolor = None, mew = 0,linestyle = '--', label = 'Model Diurnal Waveform')
    ax3.plot_date(s_time_pd.to_pydatetime(), model_s_waveform_pd, color='red',marker = 'o', markersize = 3,markeredgecolor = None, mew = 0, label = 'Model Seasonal Waveform')   
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=model_d_waveform_pd>obs_d_waveform_pd, facecolor='yellow', interpolate=True)
    ax2.fill_between(d_time_pd.to_pydatetime(),obs_d_waveform_pd, model_d_waveform_pd,where=obs_d_waveform_pd>model_d_waveform_pd, facecolor='blue', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=model_s_waveform_pd>obs_s_waveform_pd, facecolor='yellow', interpolate=True)
    ax3.fill_between(s_time_pd.to_pydatetime(),obs_s_waveform_pd, model_s_waveform_pd,where=obs_s_waveform_pd>model_s_waveform_pd, facecolor='blue', interpolate=True)
    
    ax1.set_title('Site = %s, Country = %s, Continent = %s, Process Group = %s, Lat = %s, Lon = %s, Alt = %sm,\n Data Completeness = %s%%, Anthrome Class = %s, Raw Class = %s, Grid Index = %s,%s' %(ref,country,continent,group,lat,lon,alt,complete,a_class,r_class,lat_n,lon_n))
    
    #plt.legend(loc = 'lower right') 
    plt.tight_layout()
    ax1.grid()
    ax2.grid()
    ax3.grid()
    
    pl.pop(0).remove()
    pl2.pop(0).remove()
    pl3.pop(0).remove()
    pl4.pop(0).remove()
    
    plt.show()
    
    return

    
    
def bulge_calc(obs_refs,obs_waveforms,model_waveforms,max_ph):
    #save out data to netcdf
    root_grp_period = Dataset('bulge.nc', 'w')
    root_grp_period.description = 'Bulge info found between obs and model - Program written by Dene Bowdalo'
        
    for r in range(len(obs_refs)):
        ref = obs_refs[r]
        
        print '\n',ref
        
        ref_period = root_grp_period.createGroup('%s'%(ref.lower()))
    
        obs_waveform = obs_waveforms[r]
        model_waveform = model_waveforms[r]
        waveform_len = len(obs_waveform)
        
        #process date and time, make cuts of time for day and year respectively    
        #date_l = date.astype(int)
        #time_l = time.astype(int)
        #full_times = date_process(date_l,time_l,start_year)
        #full_times_year = full_times[:8766]
        #full_times_day = full_times[:24]
        
        #get spacing
        #spacing = full_times_day[1] - full_times_day[0]
        

        #calc residual
        residual = model_waveform - obs_waveform
        
        #make ref bulge arrays
        b_sign = []
        b_maxconcoffset = []
        b_maxconcoffsetph= []
        b_aveconcoffset = []
        b_totalconcoffset = []
        b_startphaseoffset = []
        b_endphaseoffset = []
        b_totalphaseoffset = []
        b_pcbulgearea = []
        b_pcbulgetime = []
        
        #calc summary stats for period
        #) N bulges
        #) N + bulges
        #) N - bulges
        #) Ave conc. offset (+ or -) over whole period
        #) Total conc. offset over whole period (abs)
        #) Ave + conc. offset over whole period
        #) Total + area bulge over period
        #) Ave - conc. offset over whole period
        #) Total - area bulge over period
        #) Percent + area bulge of all bulges
        #) Percent - area bulge of all bulges
        #) Percent + area bulge of time
        #) Percent - area bulge of time
        
        ph_ratio = max_ph/waveform_len
        
        a_aveconcoffset = np.average(residual)
        a_totalconcoffset = np.sum(np.abs(residual))
        
        test = residual > 0
        residual_plus = residual[test]
        residual_plus_len = len(residual_plus)
        a_aveplusconcoffset = np.average(residual_plus)
        a_totalplusconcoffset = np.sum(residual_plus)
        
        test = residual < 0
        residual_neg = residual[test]
        residual_neg_len = len(residual_neg)
        a_avenegconcoffset = np.average(residual_neg)
        a_totalnegconcoffset = np.sum(residual_neg)  
        
        a_pcplusbulgearea = (100./a_totalconcoffset)*a_totalplusconcoffset
        a_pcplusbulgetime = ph_ratio * residual_plus_len
        a_pcnegbulgearea = (100./a_totalconcoffset)*np.abs(a_totalnegconcoffset)
        a_pcnegbulgetime = ph_ratio *residual_neg_len
        
        
        #find first point in residual where there is a zero crossover and iterative inds from that
        zero_crossings = np.where(np.diff(np.signbit(residual)))[0]
        
        #zero crossings will be zero if there is only 1 bulge, thus test for this
        if len(zero_crossings) > 0:
            zero_crossing_first = zero_crossings[0]+1
        else:
            zero_crossing_first = 0
        
        iter_inds = range(zero_crossing_first,len(residual),1) + range(0,zero_crossing_first,1)
        
        #calc different bulges in turn and write out stats.
        bulge_out = False
        current_sign = ''
        nbulges = 0
        nplusbulges = 0
        nnegbulges = 0
        
        last_ind = iter_inds[-1]
        
        for i in iter_inds:
            
            #if ind is last ind then finish by writing out bulge on currently
            if i == last_ind:
                save_start_i = np.copy(start_i)
                end_i = i+1
                bulge_out = True
                nbulges+=1
                if current_sign == '+':
                    nplusbulges+=1
                    bulge_flag = '+'
                else:
                    nnegbulges+=1
                    bulge_flag = '-'
                current_sign = 'na'
                print nbulges
            
            if current_sign == '+':
                if residual[i] > 0:
                    pass
                else:
                    save_start_i = np.copy(start_i)
                    end_i = i
                    bulge_out = True
                    nbulges+=1
                    nplusbulges+=1
                    current_sign = ''
                    bulge_flag = '+'
                    print nbulges
                    
            elif current_sign == '-':
                if residual[i] < 0:
                    pass
                else:
                    save_start_i = np.copy(start_i)
                    end_i = i
                    bulge_out = True
                    nbulges+=1
                    nnegbulges+=1
                    current_sign = ''
                    bulge_flag = '-'
                    print nbulges
            
            if current_sign == '':
                #get sign for next i 
                try:
                    if residual[i] > 0:
                        current_sign = '+'
                        start_i = i
                    elif residual[i] < 0:
                        current_sign = '-'
                        start_i = i
                    else:
                        current_sign = ''
                    
                except:
                    pass
                    
            if bulge_out == True:
                print bulge_flag
            
                #calc specific bulge stats
                #) Max conc offset
                #) Point of max conc offset (phase)
                #) Average conc offset
                #) Total conc offset
                #) Start phase offset
                #) End phase offset
                #) Total phase offset
                #) Percent of total bulge taken up by current bulge
                #) Percent of time taken up by bulge
        
                print save_start_i,end_i
                #make array of inds to take from residual from start_i and end_i
                if end_i < save_start_i:
                    ind_array = range(save_start_i,len(residual),1) + range(0,end_i,1)
                else:
                    ind_array = range(save_start_i,end_i,1) 
        
        
                b_sign.append(bulge_flag)
                b_maxconcoffset.append(np.max(np.abs(residual[ind_array])))
                ind_max = ind_array[np.argmax(np.abs(residual[ind_array]))]
                #print ph_ratio, ind_max
                b_maxconcoffsetph.append(ph_ratio * ind_max) 
                b_aveconcoffset.append(np.average(np.abs(residual[ind_array])))
                b_totalconcoffset.append(np.sum(np.abs(residual[ind_array])))
                b_startphaseoffset.append(ph_ratio*save_start_i)
                b_endphaseoffset.append(ph_ratio*end_i)
                b_totalphaseoffset.append(b_endphaseoffset[-1] - b_startphaseoffset[-1])
                b_pcbulgearea.append((100./a_totalconcoffset) * b_totalconcoffset[-1])
                b_pcbulgetime.append((max_ph/100.) * b_totalphaseoffset[-1])
                
                bulge_out = False
    

        for nb in range(nbulges):
            #write out from biggest total bulge to smallest 
            
            #sort all bulge arrays based on total bulge (biggest to smallest)
            sorted_inds = sorted(range(len(b_totalconcoffset)),key=lambda x:b_totalconcoffset[x])
            sorted_inds = sorted_inds[::-1]
            
            b_sign = np.array(b_sign)
            b_maxconcoffset = np.array(b_maxconcoffset)
            b_maxconcoffsetph = np.array(b_maxconcoffsetph)
            b_aveconcoffset = np.array(b_aveconcoffset)
            b_totalconcoffset = np.array(b_totalconcoffset)
            b_startphaseoffset = np.array(b_startphaseoffset)
            b_endphaseoffset = np.array(b_endphaseoffset)
            b_totalphaseoffset = np.array(b_totalphaseoffset)
            b_pcbulgearea = np.array(b_pcbulgearea)
            b_pcbulgetime = np.array(b_pcbulgetime)
            
            b_sign = b_sign[sorted_inds]
            b_maxconcoffset = b_maxconcoffset[sorted_inds]                
            b_maxconcoffsetph = b_maxconcoffsetph[sorted_inds]
            b_aveconcoffset = b_aveconcoffset[sorted_inds]
            b_totalconcoffset = b_totalconcoffset[sorted_inds]
            b_startphaseoffset = b_startphaseoffset[sorted_inds]
            b_endphaseoffset = b_endphaseoffset[sorted_inds]
            b_totalphaseoffset = b_totalphaseoffset[sorted_inds]
            b_pcbulgearea = b_pcbulgearea[sorted_inds]
            b_pcbulgetime = b_pcbulgetime[sorted_inds]
            
            bulge_grp = ref_period.createGroup('bulge%s'%(nb+1))

            bulge_grp.b_sign = b_sign[nb]
            bulge_grp.b_maxconcoffset = b_maxconcoffset[nb]
            bulge_grp.b_maxconcoffsetph = b_maxconcoffsetph[nb]
            bulge_grp.b_aveconcoffset = b_aveconcoffset[nb]
            bulge_grp.b_totalconcoffset = b_totalconcoffset[nb]
            bulge_grp.b_startphaseoffset = b_startphaseoffset[nb]
            bulge_grp.b_endphaseoffset = b_endphaseoffset[nb]
            bulge_grp.b_totalphaseoffset = b_totalphaseoffset[nb]
            bulge_grp.b_pcbulgearea = b_pcbulgearea[nb]
            bulge_grp.b_pcbulgetime = b_pcbulgetime[nb]
            

        #write out summary stats
        ref_period.nbulges = nbulges
        ref_period.nplusbulges = nplusbulges
        ref_period.nnegbulges = nnegbulges

        ref_period.a_aveconcoffset = a_aveconcoffset
        ref_period.a_aveplusconcoffset = a_aveplusconcoffset
        ref_period.a_totalplusconcoffset = a_totalplusconcoffset
        ref_period.a_avenegconcoffset = a_avenegconcoffset
        ref_period.a_totalnegconcoffset = a_totalnegconcoffset
        
        ref_period.a_pcplusbulgearea = a_pcplusbulgearea
        ref_period.a_pcnegbulgearea = a_pcnegbulgearea
        ref_period.a_pcplusbulgetime = a_pcplusbulgetime
        ref_period.a_pcnegbulgetime = a_pcnegbulgetime
        
    root_grp_period.close()
        
    return
    
def cylic_weighted_average(data,weights,max_ph):
    #print 'raw data = ',data
    
    #convert phase into radians
    pi2 = np.pi*2.
    ratio = pi2/max_ph
    rads = ratio*np.array(data)
    
    #print 'raw data rads = ',data
    
    x = y = 0.
    for rad, weight in zip(rads, weights):
        x += math.cos(rad) * weight
        y += math.sin(rad) * weight

    mean = math.atan2(y, x)

    #print 'processed data rad = ',mean

    #convert radians to 0 to 2pi (from -np.pi:pi)
    if mean < 0:
        mean = np.pi+ (np.pi-np.abs(mean))

    #convert radians back to phase
    ratio = max_ph/(np.pi*2.)
    mean = ratio*mean

    #print 'processed data = ', mean
    return mean
        
def yamartino_ave_std(data,max_ph):
    #remove nans
    data = data[~np.isnan(data)]
    
    if len(data) > 0:
        #convert phase into radians
        pi2 = np.pi*2.
        ratio = pi2/max_ph
        rads = ratio*np.array(data)
    
        sum_sin = 0
        sum_cos = 0
    
        ave = []
        
        for i in range(len(rads)):
            sum_sin = sum_sin + np.sin(rads[i])        
            sum_cos = sum_cos + np.cos(rads[i])
    
        ave = np.arctan2(sum_sin,sum_cos)
    
        #convert radians to 0 to 2pi (from -np.pi:pi)
        if ave < 0:
            ave = np.pi+ (np.pi-np.abs(ave))

        #convert radians back to phase
        ratio = max_ph/(np.pi*2.)
        mean = ratio*ave
    
        #calculate stdev
        sum_cos = sum_cos/len(data)
        sum_sin = sum_sin/len(data)
    
        e = np.sqrt(1. - sum_sin * sum_sin - sum_cos * sum_cos)
        std = np.arcsin(e) * (1 + 0.1547 * e * e * e)
        std = max_ph/(pi2/std)

        return mean,std
    else:
        return np.NaN,np.NaN
    
def maxmin_calc(d_waveforms,s_waveforms):
    #save out data to netcdf
    root_grp_period = Dataset('maxmin.nc', 'w')
    root_grp_period.description = 'MaxMin info - Program written by Dene Bowdalo'
    
    root_grp_period.createDimension('waveform_len',len(d_waveforms))
    
    d_n_peaks = []
    d_peak1_conc = []
    d_peak2_conc = []
    d_peak3_conc = []
    d_peak4_conc = []
    d_peak1_ph = []
    d_peak2_ph = []
    d_peak3_ph = []
    d_peak4_ph = []
    s_n_peaks = []
    s_peak1_conc = []
    s_peak2_conc = []
    s_peak3_conc = []
    s_peak4_conc = []
    s_peak1_ph = []
    s_peak2_ph = []
    s_peak3_ph = []
    s_peak4_ph = []
    
    print 'daily processing'
    #calculate diurnal
    max_ph = 24.
    ph_ratio = max_ph/len(d_waveforms[0])

    for wf in d_waveforms:
        peak_conc = []
        peak_ph = []
        
        n_peaks = 1
    
        max_ind = np.argmax(wf)
        iter_inds_prev = range(max_ind,len(wf),1) + range(0,max_ind-1,1)
        iter_inds_next = range(max_ind+2,len(wf),1) + range(0,max_ind+1,1)
        iter_inds = range(max_ind+1,len(wf),1) + range(0,max_ind,1)
        peak_conc.append(wf[max_ind])
        peak_ph.append(ph_ratio*max_ind)
        for i in range(len(iter_inds)):
            if (wf[iter_inds_prev[i]] < wf[iter_inds[i]]) & (wf[iter_inds_next[i]] < wf[iter_inds[i]]):
                peak_conc.append(wf[iter_inds[i]])
                peak_ph.append(ph_ratio*iter_inds[i])
                n_peaks+=1
        
        if n_peaks == 1:
            d_n_peaks.append(1)
            d_peak1_ph.append(peak_ph[0])
            d_peak2_ph.append(-99999)
            d_peak3_ph.append(-99999)
            d_peak4_ph.append(-99999)
            d_peak1_conc.append(peak_conc[0])
            d_peak2_conc.append(-99999)
            d_peak3_conc.append(-99999)
            d_peak4_conc.append(-99999)
        
        elif n_peaks == 2:
            d_n_peaks.append(2)
            d_peak1_ph.append(peak_ph[0])
            d_peak2_ph.append(peak_ph[1])
            d_peak3_ph.append(-99999)
            d_peak4_ph.append(-99999)
            d_peak1_conc.append(peak_conc[0])
            d_peak2_conc.append(peak_conc[1])
            d_peak3_conc.append(-99999)
            d_peak4_conc.append(-99999)
        
        elif n_peaks == 3:
            d_n_peaks.append(3)
            d_peak1_ph.append(peak_ph[0])
            d_peak2_ph.append(peak_ph[1])
            d_peak3_ph.append(peak_ph[2])
            d_peak4_ph.append(-99999)
            d_peak1_conc.append(peak_conc[0])
            d_peak2_conc.append(peak_conc[1])
            d_peak3_conc.append(peak_conc[2])
            d_peak4_conc.append(-99999)
        
        elif n_peaks == 4:
            d_n_peaks.append(4)
            d_peak1_ph.append(peak_ph[0])
            d_peak2_ph.append(peak_ph[1])
            d_peak3_ph.append(peak_ph[2])
            d_peak4_ph.append(peak_ph[3])
            d_peak1_conc.append(peak_conc[0])
            d_peak2_conc.append(peak_conc[1])
            d_peak3_conc.append(peak_conc[2])
            d_peak4_conc.append(peak_conc[3])
        
    print 'seasonal processing'
    #calculate seasonal
    max_ph = 12.
    ph_ratio = max_ph/len(s_waveforms[0])
    
    for wf in s_waveforms:
        n_peaks = 1
        peak_conc = []
        peak_ph = []
    
        max_ind = np.argmax(wf)
        iter_inds_prev = range(max_ind,len(wf),1) + range(0,max_ind-1,1)
        iter_inds_next = range(max_ind+2,len(wf),1) + range(0,max_ind+1,1)
        iter_inds = range(max_ind+1,len(wf),1) + range(0,max_ind,1)
        peak_conc.append(wf[max_ind])
        peak_ph.append(ph_ratio*max_ind)
        for i in range(len(iter_inds)):
            if (wf[iter_inds_prev[i]] < wf[iter_inds[i]]) & (wf[iter_inds_next[i]] < wf[iter_inds[i]]):
                peak_conc.append(wf[iter_inds[i]])
                peak_ph.append(ph_ratio*iter_inds[i])
                n_peaks+=1
        
            
        if n_peaks == 1:
            s_n_peaks.append(1)
            s_peak1_ph.append(peak_ph[0])
            s_peak2_ph.append(-99999)
            s_peak3_ph.append(-99999)
            s_peak4_ph.append(-99999)
            s_peak1_conc.append(peak_conc[0])
            s_peak2_conc.append(-99999)
            s_peak3_conc.append(-99999)
            s_peak4_conc.append(-99999)
        
        elif n_peaks == 2:
            s_n_peaks.append(2)
            s_peak1_ph.append(peak_ph[0])
            s_peak2_ph.append(peak_ph[1])
            s_peak3_ph.append(-99999)
            s_peak4_ph.append(-99999)
            s_peak1_conc.append(peak_conc[0])
            s_peak2_conc.append(peak_conc[1])
            s_peak3_conc.append(-99999)
            s_peak4_conc.append(-99999)
        
        elif n_peaks == 3:
            s_n_peaks.append(3)
            s_peak1_ph.append(peak_ph[0])
            s_peak2_ph.append(peak_ph[1])
            s_peak3_ph.append(peak_ph[2])
            s_peak4_ph.append(-99999)
            s_peak1_conc.append(peak_conc[0])
            s_peak2_conc.append(peak_conc[1])
            s_peak3_conc.append(peak_conc[2])
            s_peak4_conc.append(-99999)
        
        elif n_peaks == 4:
            s_n_peaks.append(4)
            s_peak1_ph.append(peak_ph[0])
            s_peak2_ph.append(peak_ph[1])
            s_peak3_ph.append(peak_ph[2])
            s_peak4_ph.append(peak_ph[3])
            s_peak1_conc.append(peak_conc[0])
            s_peak2_conc.append(peak_conc[1])
            s_peak3_conc.append(peak_conc[2])
            s_peak4_conc.append(peak_conc[3])
            
            
    d_n_p = root_grp_period.createVariable('diurnal_n_peaks', 'f8', ('waveform_len',))
    d_p1_ph = root_grp_period.createVariable('diurnal_peak1_phase', 'f8', ('waveform_len',))
    d_p2_ph = root_grp_period.createVariable('diurnal_peak2_phase', 'f8', ('waveform_len',))
    d_p3_ph = root_grp_period.createVariable('diurnal_peak3_phase', 'f8', ('waveform_len',))
    d_p4_ph = root_grp_period.createVariable('diurnal_peak4_phase', 'f8', ('waveform_len',))
    d_p1_conc = root_grp_period.createVariable('diurnal_peak1_conc', 'f8', ('waveform_len',))
    d_p2_conc = root_grp_period.createVariable('diurnal_peak2_conc', 'f8', ('waveform_len',))
    d_p3_conc = root_grp_period.createVariable('diurnal_peak3_conc', 'f8', ('waveform_len',))
    d_p4_conc = root_grp_period.createVariable('diurnal_peak4_conc', 'f8', ('waveform_len',))           
    s_n_p = root_grp_period.createVariable('seasonal_n_peaks', 'f8', ('waveform_len',))
    s_p1_ph = root_grp_period.createVariable('seasonal_peak1_phase', 'f8', ('waveform_len',))
    s_p2_ph = root_grp_period.createVariable('seasonal_peak2_phase', 'f8', ('waveform_len',))
    s_p3_ph = root_grp_period.createVariable('seasonal_peak3_phase', 'f8', ('waveform_len',))
    s_p4_ph = root_grp_period.createVariable('seasonal_peak4_phase', 'f8', ('waveform_len',))
    s_p1_conc = root_grp_period.createVariable('seasonal_peak1_conc', 'f8', ('waveform_len',))
    s_p2_conc = root_grp_period.createVariable('seasonal_peak2_conc', 'f8', ('waveform_len',))
    s_p3_conc = root_grp_period.createVariable('seasonal_peak3_conc', 'f8', ('waveform_len',))
    s_p4_conc = root_grp_period.createVariable('seasonal_peak4_conc', 'f8', ('waveform_len',))
    
    d_n_p[:] = d_n_peaks
    d_p1_ph[:] = d_peak1_ph
    d_p2_ph[:] = d_peak2_ph
    d_p3_ph[:] = d_peak3_ph
    d_p4_ph[:] = d_peak4_ph
    d_p1_conc[:] = d_peak1_conc
    d_p2_conc[:] = d_peak2_conc
    d_p3_conc[:] = d_peak3_conc
    d_p4_conc[:] = d_peak4_conc
    s_n_p[:] = s_n_peaks
    s_p1_ph[:] = s_peak1_ph
    s_p2_ph[:] = s_peak2_ph
    s_p3_ph[:] = s_peak3_ph
    s_p4_ph[:] = s_peak4_ph
    s_p1_conc[:] = s_peak1_conc
    s_p2_conc[:] = s_peak2_conc
    s_p3_conc[:] = s_peak3_conc
    s_p4_conc[:] = s_peak4_conc

        
    root_grp_period.close()
        
    return
  
def area_dicts():  
    area_boundaries = {'CE_NA':[37,45,-90,-50],'SE_NA':[25,37,-90,-60],'S_NA':[25,37,-105,-90],'SW_NA':[25,40,-130,-105],'C_NA':[37,70,-115,-90],'NE_NA':[45,80,-90,-50],'NW_NA':[40,80,-141,-115],'AL':[50,80,-170,-141],'N_EU':[55,80,5,45],'C_EU':[47,55,4,16],'CS_EU':[44,47,4,16],'CW_EU':[44,50,-6,4],'SW_EU':[35,44,-10,-1],'S_EU':[30,44,-1,40],'E_EU':[44,55,16,40],'NW_EU':[50,70,-15,2],'NE_AS':[30,60,90,155],'SE_AS':[-10,30,90,130],'C_AS':[30,60,60,90],'S_AS':[0,30,60,90],'N_O':[0,0,0,0],'S_O':[0,0,0,0],'OC':[0,0,0,0],'AF':[0,0,0,0],'SA':[0,0,0,0],'ANT':[0,0,0,0],'ARC':[0,0,0,0]}
    area_tag = {'CE_NA':['NA'],'SE_NA':['NA'],'S_NA':['NA'],'SW_NA':['NA'],'C_NA':['NA'],'NE_NA':['NA'],'NW_NA':['NA'],'AL':['NA'],'C_EU':['EU'],'N_EU':['EU'],'S_EU':['EU'],'E_EU':['EU'],'NW_EU':['EU'],'NE_AS':['AS'],'SE_AS':['AS'],'C_AS':['AS'],'S_AS':['AS'],'N_O':['O'],'S_O':['O'],'OC':['OC'],'AF':['AF'],'SA':['SA'],'ANT':['ANT'],'ARC':['ARC']}
    area_labels = {'CE_NA':'CE NA','SE_NA':'SE NA','S_NA':'S NA','SW_NA':'SW NA','C_NA':'C NA','NE_NA':'NE NA','NW_NA':'NW NA','AL':'Alaska','N_EU':'N EU','NW_EU':'NW EU','C_EU':'C EU','E_EU':'E EU','S_EU':'S EU','NE_AS':'NE Asia','SE_AS':'SE Asia','C_AS':'C Asia','S_AS':'S Asia','N_O':'NH Oceanic','S_O':'SH Oceanic','OC':'Oceania','AF':'Africa','SA':'S America','ANT':'Antarctica','ARC':'Arctic'}
    
    return area_boundaries,area_tag,area_labels
    
def area_cut(area,obs_lats,obs_lons,tags,area_grid,area_tag):
    if (area == 'C_EU') or (area == 'E_EU') or (area == 'S_EU') or (area == 'NW_EU') or (area == 'N_EU') or (area == 'CW_EU') or (area == 'SW_EU') or (area == 'CS_EU'):
        cut_test = (obs_lats >= area_grid[0]) & (obs_lats < area_grid[1]) & (obs_lons >= area_grid[2]) & (obs_lons < area_grid[3]) & (tags == area_tag[0])
    elif (area == 'OC') or (area == 'AF') or (area == 'ANT') or (area == 'ARC') or (area == 'SA'):
        cut_test = tags == area_tag[0]
    elif (area == 'N_O'):
        cut_test = (tags == area_tag[0]) & (obs_lats >= 0)
    elif (area == 'S_O'):
        cut_test = (tags == area_tag[0]) & (obs_lats < 0)
    elif (area == 'CE_NA') or (area == 'SE_NA') or (area == 'S_NA') or (area == 'SW_NA') or (area == 'C_NA') or (area == 'NE_NA') or (area == 'NW_NA') or (area == 'AL'):              
        cut_test = (obs_lats >= area_grid[0]) & (obs_lats < area_grid[1]) & (obs_lons >= area_grid[2]) & (obs_lons < area_grid[3]) & (tags == area_tag[0])
    elif (area == 'NE_AS') or (area == 'SE_AS') or (area == 'S_AS') or (area == 'C_AS'):
        cut_test = (obs_lats >= area_grid[0]) & (obs_lats < area_grid[1]) & (obs_lons >= area_grid[2]) & (obs_lons < area_grid[3]) & (tags == area_tag[0])

    return cut_test
    
def get_area(area_boundaries,tag,obs_lat,obs_lon,obs_ref):
    
    if tag == 'NA':
        if (obs_lat >= area_boundaries['CE_NA'][0]) & (obs_lat < area_boundaries['CE_NA'][1]) & (obs_lon >= area_boundaries['CE_NA'][2]) & (obs_lon < area_boundaries['CE_NA'][3]):
            area = 'CE_NA'
        elif (obs_lat >= area_boundaries['SE_NA'][0]) & (obs_lat < area_boundaries['SE_NA'][1]) & (obs_lon >= area_boundaries['SE_NA'][2]) & (obs_lon < area_boundaries['SE_NA'][3]):
            area = 'SE_NA'
        elif (obs_lat >= area_boundaries['S_NA'][0]) & (obs_lat < area_boundaries['S_NA'][1]) & (obs_lon >= area_boundaries['S_NA'][2]) & (obs_lon < area_boundaries['S_NA'][3]):
            area = 'S_NA'
        elif (obs_lat >= area_boundaries['SW_NA'][0]) & (obs_lat < area_boundaries['SW_NA'][1]) & (obs_lon >= area_boundaries['SW_NA'][2]) & (obs_lon < area_boundaries['SW_NA'][3]):
            area = 'SW_NA'
        elif (obs_lat >= area_boundaries['C_NA'][0]) & (obs_lat < area_boundaries['C_NA'][1]) & (obs_lon >= area_boundaries['C_NA'][2]) & (obs_lon < area_boundaries['C_NA'][3]):
            area = 'C_NA'
        elif (obs_lat >= area_boundaries['NE_NA'][0]) & (obs_lat < area_boundaries['NE_NA'][1]) & (obs_lon >= area_boundaries['NE_NA'][2]) & (obs_lon < area_boundaries['NE_NA'][3]):
            area = 'NE_NA'
        elif (obs_lat >= area_boundaries['NW_NA'][0]) & (obs_lat < area_boundaries['NW_NA'][1]) & (obs_lon >= area_boundaries['NW_NA'][2]) & (obs_lon < area_boundaries['NW_NA'][3]):
            area = 'NW_NA'
        elif (obs_lat >= area_boundaries['AL'][0]) & (obs_lat < area_boundaries['AL'][1]) & (obs_lon >= area_boundaries['AL'][2]) & (obs_lon < area_boundaries['AL'][3]):
            area = 'AL'
        else:
            area = 'NO GIVEN AREA'
            print tag,obs_lat,obs_lon,obs_ref
            
    elif tag == 'EU':
        if (obs_lat >= area_boundaries['N_EU'][0]) & (obs_lat < area_boundaries['N_EU'][1]) & (obs_lon >= area_boundaries['N_EU'][2]) & (obs_lon < area_boundaries['N_EU'][3]):
            area = 'N_EU'
        elif (obs_lat >= area_boundaries['C_EU'][0]) & (obs_lat < area_boundaries['C_EU'][1]) & (obs_lon >= area_boundaries['C_EU'][2]) & (obs_lon < area_boundaries['C_EU'][3]):
            area = 'C_EU'
        elif (obs_lat >= area_boundaries['CS_EU'][0]) & (obs_lat < area_boundaries['CS_EU'][1]) & (obs_lon >= area_boundaries['CS_EU'][2]) & (obs_lon < area_boundaries['CS_EU'][3]):
            area = 'CS_EU'  
        elif (obs_lat >= area_boundaries['CW_EU'][0]) & (obs_lat < area_boundaries['CW_EU'][1]) & (obs_lon >= area_boundaries['CW_EU'][2]) & (obs_lon < area_boundaries['CW_EU'][3]):
            area = 'CW_EU'
        elif (obs_lat >= area_boundaries['SW_EU'][0]) & (obs_lat < area_boundaries['SW_EU'][1]) & (obs_lon >= area_boundaries['SW_EU'][2]) & (obs_lon < area_boundaries['SW_EU'][3]):
            area = 'SW_EU'
        elif (obs_lat >= area_boundaries['S_EU'][0]) & (obs_lat < area_boundaries['S_EU'][1]) & (obs_lon >= area_boundaries['S_EU'][2]) & (obs_lon < area_boundaries['S_EU'][3]):
            area = 'S_EU'
        elif (obs_lat >= area_boundaries['NW_EU'][0]) & (obs_lat < area_boundaries['NW_EU'][1]) & (obs_lon >= area_boundaries['NW_EU'][2]) & (obs_lon < area_boundaries['NW_EU'][3]):
            area = 'NW_EU'
        elif (obs_lat >= area_boundaries['E_EU'][0]) & (obs_lat < area_boundaries['E_EU'][1]) & (obs_lon >= area_boundaries['E_EU'][2]) & (obs_lon < area_boundaries['E_EU'][3]):
            area = 'E_EU'
        else:
            area = 'NO GIVEN AREA'
            print tag,obs_lat,obs_lon,obs_ref
     
    elif tag == 'AS':
        if (obs_lat >= area_boundaries['NE_AS'][0]) & (obs_lat < area_boundaries['NE_AS'][1]) & (obs_lon >= area_boundaries['NE_AS'][2]) & (obs_lon < area_boundaries['NE_AS'][3]):
            area = 'NE_AS'
        elif (obs_lat >= area_boundaries['SE_AS'][0]) & (obs_lat < area_boundaries['SE_AS'][1]) & (obs_lon >= area_boundaries['SE_AS'][2]) & (obs_lon < area_boundaries['SE_AS'][3]):
            area = 'SE_AS'    
        elif (obs_lat >= area_boundaries['C_AS'][0]) & (obs_lat < area_boundaries['C_AS'][1]) & (obs_lon >= area_boundaries['C_AS'][2]) & (obs_lon < area_boundaries['C_AS'][3]):
            area = 'C_AS'
        elif (obs_lat >= area_boundaries['S_AS'][0]) & (obs_lat < area_boundaries['S_AS'][1]) & (obs_lon >= area_boundaries['S_AS'][2]) & (obs_lon < area_boundaries['S_AS'][3]):
            area = 'S_AS'
        else:
            area = 'NO GIVEN AREA'
            print tag,obs_lat,obs_lon,obs_ref
        
    elif tag == 'OC':
        area = 'OC'
        
    elif tag == 'ANT':
        area = 'ANT'
      
    elif tag == 'ARC':
        area = 'ARC'  
    
    elif tag == 'AF':
        area = 'AF'
    
    elif tag == 'SA':
        area = 'SA'
    
    elif (tag == 'O') & (obs_lat >= 0):
        area = 'N_O'
        
    elif (tag == 'O') & (obs_lat < 0):
        area = 'S_O'
                
    return area
    
def get_pressure_edges(psfc):

    # !  GEOS-4, GEOS-5, GEOS-5.7, and MERRA (hybrid grids):
    # !  ----------------------------------------------------------------------------
    # !  For GEOS-4/GEOS-5/MERRA met data products, the pressure at the bottom edge 
    # !  of grid box (I,J,L) is defined as follows:
    # !                                                                             .
    # !     Pedge(I,J,L) = Ap(L) + [ Bp(L) * Psurface(I,J) ]
    # !                                                                             .
    # !  where
    # !                                                                             .
    # !     Psurface(I,J) is  the "true" surface pressure at lon,lat (I,J)
    # !     Ap(L)         has the same units as surface pressure [hPa]
    # !     Bp(L)         is  a unitless constant given at level edges
    # !                                                                             .
    # !  Ap(L) and Bp(L) are given to us by GMAO.
    # 
    # 
    # !  The following are true for GCAP, GEOS-3, GEOS-4, GEOS-5, MERRA:
    # !  ----------------------------------------------------------------------------
    # !  (1) Bp(LLPAR+1) = 0.0          (L=LLPAR+1 is the atmosphere top)
    # !  (2) Bp(1)       = 1.0          (L=1       is the surface       )
    # !  (3) PTOP        = Ap(LLPAR+1)  (L=LLPAR+1 is the atmosphere top) 
    # 
    # 
    # !-----------------------------------------------------------------
    # ! GEOS-5/MERRA 47-level reduced vertical grid
    # !  
    # !  Bottom   Bottom    # levels
    # !  edge of  edge prs  lumped 
    # !  level    (hPa)     together
    # !
    # !   PTOP       0.010   
    # !    47        0.066     4
    # !    46        0.211     4
    # !    45        0.617     4
    # !    44        1.651     4
    # !    43        4.077     4
    # !    42        9.293     4
    # !    41       19.792     4
    # !    40       28.368     2
    # !    39       40.175     2
    # !    38       56.388     2
    # !    37       78.512     2
    # ! %%%% START LUMPING LEVELS ABOVE HERE %%%%%
    # !    36       92.366       
    # !    35      108.663
    # !    34      127.837
    # !    33      150.393
    # !    32      176.930
    # ! %%%% FIXED-PRESSURE LEVELS BEGIN HERE %%%%
    # !-----------------------------------------------------------------

    # Ap [hPa] for 47 levels (48 edges)
    AP = np.array([0.000000e+00, 4.804826e-02, 6.593752e+00, 1.313480e+01,
            1.961311e+01, 2.609201e+01, 3.257081e+01, 3.898201e+01,
            4.533901e+01, 5.169611e+01, 5.805321e+01, 6.436264e+01,
            7.062198e+01, 7.883422e+01, 8.909992e+01, 9.936521e+01,
            1.091817e+02, 1.189586e+02, 1.286959e+02, 1.429100e+02,
            1.562600e+02, 1.696090e+02, 1.816190e+02, 1.930970e+02,
            2.032590e+02, 2.121500e+02, 2.187760e+02, 2.238980e+02,
            2.243630e+02, 2.168650e+02, 2.011920e+02, 1.769300e+02,
            1.503930e+02, 1.278370e+02, 1.086630e+02, 9.236572e+01,
            7.851231e+01, 5.638791e+01, 4.017541e+01, 2.836781e+01, 
            1.979160e+01, 9.292942e+00, 4.076571e+00, 1.650790e+00, 
            6.167791e-01, 2.113490e-01, 6.600001e-02, 1.000000e-02])

    # Bp [unitless] for 47 levels (48 edges)
    BP = np.array([1.000000e+00, 9.849520e-01, 9.634060e-01, 9.418650e-01,
            9.203870e-01, 8.989080e-01, 8.774290e-01, 8.560180e-01,
            8.346609e-01, 8.133039e-01, 7.919469e-01, 7.706375e-01,
            7.493782e-01, 7.211660e-01, 6.858999e-01, 6.506349e-01,
            6.158184e-01, 5.810415e-01, 5.463042e-01, 4.945902e-01,
            4.437402e-01, 3.928911e-01, 3.433811e-01, 2.944031e-01,
            2.467411e-01, 2.003501e-01, 1.562241e-01, 1.136021e-01,
            6.372006e-02, 2.801004e-02, 6.960025e-03, 8.175413e-09,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00,
            0.000000e+00, 0.000000e+00, 0.000000e+00, 0.000000e+00])

    pressure_edges = np.empty((psfc.shape[0],48,psfc.shape[1],psfc.shape[2]))
    lat_inds = range(psfc.shape[1])
    lon_inds = range(psfc.shape[2])
    for i in range(psfc.shape[0]):
        print i
        for lat_i in lat_inds:
            for lon_i in lon_inds:
                pressure_edges[i,:,lat_i,lon_i] = AP + (BP * psfc[i,lat_i,lon_i])
        print psfc[i,lat_i,lon_i]
        print pressure_edges[i,:,lat_i,lon_i]
    
    return pressure_edges
    
def quality_check(full_data,data_valid,data_resolution,alt,grid_dates,start_year,end_year):

    #check that data resolution 
    #do by checking checking differences with 1st valid point and all other valid points
    #exclude 0 differences, if minimum difference is > data_resolution min exclude site as data resolution is too poor
    test = full_data != -99999 
    val_data = full_data[test]  
    min_diff = np.abs(np.diff(val_data))
    test = min_diff != 0
    min_diff = min_diff[test]
    try:
        min_diff = np.min(min_diff)    
        if min_diff > data_resolution:
            data_valid = False
            print 'Data resolution is not sufficient.'
            print 'min diff = ', min_diff
    except:
        data_valid = False
        print 'Data resolution is not sufficient.'

    #-------------------------------------------------------    
    #Remove sites above 1000m from sea level
    if np.float64(alt) >= 1000:
        data_valid = False
        print 'Site is over 1000m from sea level.'    
        
    #-------------------------------------------------------
    #check there are no continuous gaps > 2 months, in each year which has data
    year_inv_count = 0
    year_range = range(start_year,end_year+1)
    all_years = [i[:4] for i in grid_dates]
    all_years = np.array(all_years)
    for i in range(len(year_range)):
        year_test = all_years == str(year_range[i])
        year_vals = full_data[year_test]
        
        #if all values in year are empty then continue to next year
        all_invalid = all(x == year_vals[0] for x in year_vals)
        if all_invalid == True:
            continue
        
        max_count = 0
        count = 0
        for j in year_vals:
            if j == -99999:
                count+=1
            else:   
                count = 0
            
            if count > max_count:
                max_count = count
        
        # if there is a gap > 60 days, make whole year invalid
        if max_count > 1440:
            full_data[year_test] = -99999
    
    return full_data,data_valid

def quality_check_nr(full_data,data_valid,data_resolution,alt,grid_dates,start_year,end_year):

    #check that data resolution 
    #do by checking checking differences with 1st valid point and all other valid points
    #exclude 0 differences, if minimum difference is > data_resolution min exclude site as data resolution is too poor
    test = full_data != -99999 
    val_data = full_data[test]  
    min_diff = np.abs(np.diff(val_data))
    test = min_diff != 0
    min_diff = min_diff[test]
    try:
        min_diff = np.min(min_diff)    
        if min_diff > data_resolution:
            data_valid = False
            print 'Data resolution is not sufficient.'
            print 'min diff = ', min_diff
    except:
        data_valid = False
        print 'Data resolution is not sufficient.'
        
    #-------------------------------------------------------
    #check there are no continuous gaps > 2 months, in each year which has data
    year_inv_count = 0
    year_range = range(start_year,end_year+1)
    all_years = [i[:4] for i in grid_dates]
    all_years = np.array(all_years)
    for i in range(len(year_range)):
        year_test = all_years == str(year_range[i])
        year_vals = full_data[year_test]
        
        #if all values in year are empty then continue to next year
        all_invalid = all(x == year_vals[0] for x in year_vals)
        if all_invalid == True:
            continue
        
        max_count = 0
        count = 0
        for j in year_vals:
            if j == -99999:
                count+=1
            else:   
                count = 0
            
            if count > max_count:
                max_count = count
        
        # if there is a gap > 60 days, make whole year invalid
        if max_count > 1440:
            full_data[year_test] = -99999
    
    return full_data,data_valid
    
def quality_check_periodic(full_data,data_valid,data_resolution,alt,grid_dates,start_year,end_year):

    #check that data resolution 
    #do by checking checking differences with 1st valid point and all other valid points
    #exclude 0 differences, if minimum difference is > data_resolution min exclude site as data resolution is too poor
    test = full_data != -99999 
    val_data = full_data[test]  
    min_diff = np.abs(np.diff(val_data))
    test = min_diff != 0
    min_diff = min_diff[test]
    try:
        min_diff = np.min(min_diff)    
        if min_diff > data_resolution:
            data_valid = False
            print 'Data resolution is not sufficient.'
            print 'min diff = ', min_diff
    except:
        data_valid = False
        print 'Data resolution is not sufficient.'

    #-------------------------------------------------------    
    #Remove sites above 1000m from sea level
    if np.float64(alt) >= 1000:
        data_valid = False
        print 'Site is over 1000m from sea level.'    
    
    #-------------------------------------------------------
    #check amount of valid data
    a_test = full_data != -99999
    valid = full_data[a_test]
    print len(full_data)
    print len(valid)
    valid_line = len(full_data) /2
    data_complete = (100./len(full_data)) * len(valid)
    if len(valid) < valid_line:
        data_valid = False
        print 'Site Invalid, only %s %% of data valid'%(data_complete)
    
    #-------------------------------------------------------
    #check there are no gaps > 365 days
    count = 0
    max_count = 0
    for i in full_data:
        if i == -99999:
            count+=1
        else:
            count = 0

        if count > max_count:
            max_count = count

    if max_count > 8759:
        data_valid = False
        print 'Site Invalid, gaps greater than 365 days in data.'
    
    #-------------------------------------------------------
    #check there are no continuous gaps > 2 months, in each year which has data
    year_inv_count = 0
    year_range = range(start_year,end_year+1)
    all_years = [i[:4] for i in grid_dates]
    all_years = np.array(all_years)
    for i in range(len(year_range)):
        year_test = all_years == str(year_range[i])
        year_vals = full_data[year_test]
        
        max_count = 0
        count = 0
        for j in year_vals:
            if j == -99999:
                count+=1
            else:   
                count = 0
            
            if count > max_count:
                max_count = count
        
        # if there is a gap > 60 days, year is invalid
        if max_count > 1440:
            year_inv_count+=1
    
    if year_inv_count >= int(round(len(all_years)/2.)):
        data_valid = False
        print 'Persisent Data gap in at least half of years > 60 days'
    
    return full_data,data_valid,data_complete
    
def take_average(full_data,obs_time_pd,output_res):

    #make sure all invalids set to np.nan
    invalid_test = full_data < 0
    full_data[invalid_test] = np.nan
    
    #take average of dataset
    obs_var_pd = pd.Series(full_data, index=obs_time_pd)

    if output_res == 'D':
        obs_key = lambda x: pd.Period(str(x.year)+'-'+str(x.month)+'-'+str(x.day))
    elif output_res == 'M':
        obs_key = lambda x: pd.Period(str(x.year)+'-'+str(x.month))

    #group the data
    group=obs_var_pd.groupby(obs_key)
    datafrac=group.count().div(group.size(),axis='rows')

    #calculate the mean
    obs_mean = group.mean()
    #now, zero out the mean data where datafrac<0.2
    obs_mean[datafrac<0.2]=np.NaN
    full_data = np.array(obs_mean.values.tolist())

    invalid_test = np.isnan(full_data)
    full_data[invalid_test] = -99999
    
    #------------------------------------------------------
    #check array has some valid values otherwise do not write out
    all_invalid = all(x == full_data[0] for x in full_data)
    if all_invalid == True:
        has_data = False
    else:
        has_data = True
        
    return full_data, has_data
    
def write_out_data(site_ref,process_group,root_grp,species,full_data,output_res,lat,lon,alt,raw_class_name,anthrome_class_name,mm,st,file_res,output_res_times,obs_time_pd,data_valid):
    if data_valid == True: 
        country = ''
        while country == '':
            try:
                country = get_country(lat,lon)                       
            except:
                pass

        if (output_res == 'D') or (output_res == 'M'):
            #average data
            full_data,has_data = take_average(full_data,obs_time_pd,output_res)
        
        elif (output_res == 'H'):
            #check site has valid data
            all_invalid = all(x == full_data[0] for x in full_data)
            if all_invalid == True:
                has_data = False
            else:
                has_data = True

        if has_data == True:
            #save out netcdf file
            ref = root_grp.createGroup('%s'%(site_ref.lower()))

            #set variables
            times = ref.createVariable('time', 'f8', ('time',))
            spec = ref.createVariable(species.lower(), 'f4', ('species',))

            #set group attributes
            ref.latitude = lat
            ref.longitude = lon
            ref.altitude = alt
            ref.process_group = process_group
            ref.country = country
            ref.anthrome_site_class = anthrome_class_name
            ref.raw_site_class = raw_class_name
            ref.unit = 'ppbv'
            ref.native_resolution =  file_res
            ref.measurement_method = mm
            ref.sampling_type = st

            times[:] = output_res_times
            spec[:] = full_data
            print 'Site is Valid\n'
        else:
            print 'Not writing as all monthly values are invalid'
            print 'Site is Not Valid\n'
    else:
        print 'Site is Not Valid\n' 
        
        
def write_out_data_periodic(site_ref,process_group,root_grp,species,full_data,output_res,lat,lon,alt,raw_class_name,anthrome_class_name,mm,st,file_res,data_complete,output_res_times,obs_time_pd,data_valid):
    if data_valid == True: 
        country = ''
        while country == '':
            try:
                country = get_country(lat,lon)                       
            except:
                pass

        if (output_res == 'D') or (output_res == 'M'):
            #average data
            full_data,has_data = take_average(full_data,obs_time_pd,output_res)
        
        elif (output_res == 'H'):
            #check site has valid data
            all_invalid = all(x == full_data[0] for x in full_data)
            if all_invalid == True:
                has_data = False
            else:
                has_data = True

        if has_data == True:
            #save out netcdf file
            ref = root_grp.createGroup('%s'%(site_ref.lower()))

            #set variables
            times = ref.createVariable('time', 'f8', ('time',))
            spec = ref.createVariable(species.lower(), 'f4', ('species',))

            #set group attributes
            ref.latitude = lat
            ref.longitude = lon
            ref.altitude = alt
            ref.process_group = process_group
            ref.country = country
            ref.anthrome_site_class = anthrome_class_name
            ref.raw_site_class = raw_class_name
            ref.unit = 'ppbv'
            ref.native_resolution =  file_res
            ref.measurement_method = mm
            ref.sampling_type = st
            ref.data_completeness = data_complete
            
            print len(output_res_times)
            print output_res_times
            print times
            times[:] = output_res_times
            spec[:] = full_data
            print 'Site is Valid\n'
        else:
            print 'Not writing as all monthly values are invalid'
            print 'Site is Not Valid\n'
    else:
        print 'Site is Not Valid\n' 

def read_obs_one(obs_file,species,ref,start_year,end_year):
    root_grp = Dataset(obs_file)   
    
    ref_group = root_grp.groups[ref]
    std_var = ref_group.variables[species.lower()][:]
    obs_lat = ref_group.latitude
    obs_lon = ref_group.longitude
    obs_alt = ref_group.altitude
    raw_time = ref_group.variables['time'][:]
    
    #cut time
    start_time = date2num(datetime.datetime(int(start_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    end_time = date2num(datetime.datetime(int(end_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    start_ind = np.searchsorted(raw_time,start_time)
    end_ind = np.searchsorted(raw_time,end_time)
    if (start_time < raw_time[0]) & (end_time > raw_time[-1]):
        raw_time = raw_time[:]
        std_var = std_var[:]
    elif start_time < raw_time[0]:
        raw_time = raw_time[:end_ind]
        std_var = std_var[:end_ind]
    elif end_time > raw_time[-1]:
        raw_time = raw_time[start_ind:]
        std_var = std_var[start_ind:]
    else:
        raw_time = raw_time[start_ind:end_ind]
        std_var = std_var[start_ind:end_ind]
    
    ref_time = raw_time - raw_time[0]
    
    datetime_time = num2date(raw_time,units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    
    #get raw and ref time in days
    raw_time = raw_time / 24.
    ref_time = ref_time / 24.
    
    std_var_mask = np.ma.masked_where(std_var<=0,std_var)
    pd_var = pd.Series(std_var_mask, index=datetime_time)
    
    return refs,raw_time,ref_time,datetime_time,std_var,pd_var,obs_lat,obs_lon,obs_alt 
 
def read_obs_all(obs_file,species,start_year,end_year):
    
    root_grp = Dataset(obs_file)
    refs_dict = root_grp.groups
    
    refs = []
    obs_lats = []
    obs_lons = []
    obs_alt = []
    std_var = []
    gap_inds = []
    
    for i in refs_dict.keys():
        i = i.encode('ascii')
        refs = np.append(refs,i)
    
    for ref in refs:
        ref_group = root_grp.groups[ref]
        std_var.append(ref_group.variables[species.lower()][:])
        obs_lats = np.append(obs_lats,ref_group.latitude)
        obs_lons = np.append(obs_lons,ref_group.longitude)
        obs_alt = np.append(obs_alt,ref_group.altitude)
        gap_inds.append(std_var < 0)
        
    raw_time = ref_group.variables['time'][:]
    
    start_time = date2num(datetime.datetime(int(start_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    end_time = date2num(datetime.datetime(int(end_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    start_ind = np.searchsorted(raw_time,start_time)
    end_ind = np.searchsorted(raw_time,end_time)
    if (start_time < raw_time[0]) & (end_time > raw_time[-1]):
        raw_time = raw_time[:]
        std_var = std_var[:]
    elif start_time < raw_time[0]:
        raw_time = raw_time[:end_ind]
        std_var = std_var[:end_ind]
    elif end_time > raw_time[-1]:
        raw_time = raw_time[start_ind:]
        std_var = std_var[start_ind:]
    else:
        raw_time = raw_time[start_ind:end_ind]
        std_var = std_var[start_ind:end_ind]
    
    ref_time = raw_time - raw_time[0]
    datetime_time = num2date(raw_time,units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    
    #get raw and ref time in days
    raw_time = raw_time / 24.
    ref_time = ref_time / 24.
    
    for i in range(len(refs)):
        refs[i] = refs[i].lower()

    return refs,raw_time,ref_time,datetime_time,std_var,obs_lats,obs_lons,obs_alt,gap_inds
 
def read_model_one(model_file,species,obs_lat,obs_lon,start_year,end_year):

    root_grp = Dataset(model_file)
    std_var = root_grp.variables[species.lower()][:]
    raw_time = root_grp.variables['time'][:]
    lat_e = root_grp.variables['lat_edges'][:]
    lon_e = root_grp.variables['lon_edges'][:]
    lat_c = root_grp.variables['lat_centre'][:]
    lon_c = root_grp.variables['lon_centre'][:]
    grid_size = root_grp.variables['grid_size'][:]
    grid_size = grid_size[0]

    gridbox_count = len(lat_c)*len(lon_c)
    
    start_time = date2num(datetime.datetime(int(start_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    end_time = date2num(datetime.datetime(int(end_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    start_ind = np.searchsorted(raw_time,start_time)
    end_ind = np.searchsorted(raw_time,end_time)
    if (start_time < raw_time[0]) & (end_time > raw_time[-1]):
        raw_time = raw_time[:]
        std_var = std_var[:]
    elif start_time < raw_time[0]:
        raw_time = raw_time[:end_ind]
        std_var = std_var[:end_ind]
    elif end_time > raw_time[-1]:
        raw_time = raw_time[start_ind:]
        std_var = std_var[start_ind:]
    else:
        raw_time = raw_time[start_ind:end_ind]
        std_var = std_var[start_ind:end_ind]
    
    ref_time = raw_time - raw_time[0]
    datetime_time = num2date(raw_time,units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    
    #get raw and ref time in days
    raw_time = raw_time / 24.
    ref_time = ref_time / 24.
    
    lat_n,lon_n = modules.obs_model_gridbox(lat_e,lon_e,obs_lat,obs_lon)
    std_var = std_var[:,lat_n,lon_n]
    std_var_mask = np.ma.masked_where(std_var<=0,std_var)
    pd_var = pd.Series(std_var_mask, index=datetime_time)

    return raw_time,ref_time,datetime_time,std_var,pd_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count
     

def read_model_all(model_file,species,start_year,end_year):

    root_grp = Dataset(model_file)
    std_var = root_grp.variables[species.lower()][:]
    raw_time = root_grp.variables['time'][:]
    lat_e = root_grp.variables['lat_edges'][:]
    lon_e = root_grp.variables['lon_edges'][:]
    lat_c = root_grp.variables['lat_centre'][:]
    lon_c = root_grp.variables['lon_centre'][:]
    grid_size = root_grp.variables['grid_size'][:]
    grid_size = grid_size[0]

    gridbox_count = len(lat_c)*len(lon_c)
    
    start_time = date2num(datetime.datetime(int(start_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    end_time = date2num(datetime.datetime(int(end_year),1,1,0,0,0,0),units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    start_ind = np.searchsorted(raw_time,start_time)
    end_ind = np.searchsorted(raw_time,end_time)
    if (start_time < raw_time[0]) & (end_time > raw_time[-1]):
        raw_time = raw_time[:]
        std_var = std_var[:]
    elif start_time < raw_time[0]:
        raw_time = raw_time[:end_ind]
        std_var = std_var[:end_ind]
    elif end_time > raw_time[-1]:
        raw_time = raw_time[start_ind:]
        std_var = std_var[start_ind:]
    else:
        raw_time = raw_time[start_ind:end_ind]
        std_var = std_var[start_ind:end_ind]
    
    ref_time = raw_time - raw_time[0]
    datetime_time = num2date(raw_time,units='hours since 0001-01-01 00:00:00', calendar='gregorian')
    
    #get raw and ref time in days
    raw_time = raw_time / 24.
    ref_time = ref_time / 24.
    
    return raw_time,ref_time,datetime_time,std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count


    