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
from matplotlib.ticker import FormatStrFormatter
from matplotlib.ticker import LinearLocator

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



def NOAA_data_reader_mace_head(filename):
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

def savitzky_golay(y, window_size, order, deriv=0, rate=1):
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
        ofac= 2.0001

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
    
    def readfile_gaw(filename, location):

        read = np.load(filename)
        names = read[0,1:]
        locs = np.where(read[:,0] == location)
        big = read[locs]
        print big
        valid_data = big[:,1:]
        names = names.tolist()
        valid_data = np.float64(valid_data)

        return valid_data, names

    try: 
        names
    except NameError:
# Readin the model output

        model , names = readfile("GEOS_v90103_nested_europe_GAW_logs.npy","112") #112 represents Mace Head
        model2, gaw_names = readfile_gaw("gaw_logs_O3.npy","112")  

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

# Processes the gaw baseline date 
        year=(model2[:,0]//10000)
        month=((model2[:,0]-year*10000)//100)
        day=(model2[:,0]-year*10000-month*100)

        hour=model2[:,1]//100
        min=(model2[:,1]-hour*100)

        doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(2006,1,1,0,0,0) \
              for i in range(len(year))]

        since2006_gaw=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]


#now read in the observations

    date, time, vals =  NOAA_data_reader_mace_head('O3_mace_head_ppbV.txt')

# OK need to conver values from a list to a numpy array

    valid = vals >= 0
    vals = vals[valid]
    date = date[valid]
    time = time[valid]

# Processes the date 
    year=(date//10000)
    month=((date-year*10000)//100)
    day=(date-year*10000-month*100)

    hour=time//100
    min=(time-hour*100)

    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(2006,1,1,0,0,0) \
              for i in range(len(year))]

    since2006_obs=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]
    
    since2006_obs = np.array(since2006_obs) 
    since2006_obs_gaw = since2006_obs 
    

    valids2= since2006_obs <= np.max(since2006)
    since2006_obs=since2006_obs[valids2]
    var3=vals[valids2]

    valids3= since2006_obs >= np.min(since2006)
    since2006_obs=since2006_obs[valids3]
    var4=var3[valids3]

#Pre normalise obs data for lomb analysis
    standard_deviation_obs_p = np.std(var4)
    mean_obs_p = np.mean(var4)
    normal_var2 = var4-mean_obs_p
    normal_var2 = normal_var2/standard_deviation_obs_p
 
#Pre normalise obs data for lomb analysis
    standard_deviation_obs_p_gaw = np.std(vals)
    mean_obs_p_gaw = np.mean(vals)
    normal_var2_gaw = vals-mean_obs_p_gaw
    normal_var2_gaw = normal_var2_gaw/standard_deviation_obs_p_gaw

    print 'obs', normal_var2_gaw
#Calculate variance of pre-processed obs data- should be 1 if normal
    #standard_dev_obs = np.std(normal_var_2, dtype=np.float64)
    #variance_obs = standard_dev_obs**2
    #print 'Variance - pre-processed obs data= ', variance_obs

#Define sampling intervals
    samp_spacing = 1./24.

#Convert model time array into numpy array
    since2006=np.array(since2006)
    since2006_gaw=np.array(since2006_gaw)
#Need to normalise model data also
    if model_cut_switch == 0:
        k=names.index(species)
        model_cut = model[:,k]*unit_cut
        print 'model cut', model_cut  	
    if model_cut_switch == 1:
        model_cut = read_diff_species()
    standard_deviation_model_p = np.std(model_cut)
    mean_model_p = np.mean(model_cut)
    normal_model = model_cut-mean_model_p
    normal_model = normal_model/standard_deviation_model_p
    print 'normal model', normal_model

#Need to normalise model baseline data also
    j=gaw_names.index(species)
    model_cut_gaw = model2[:,j]*unit_cut
    standard_deviation_model_p_gaw = np.std(model_cut_gaw)
    mean_model_p_gaw = np.mean(model_cut_gaw)
    normal_model_gaw = model_cut_gaw-mean_model_p_gaw
    normal_model_gaw = normal_model_gaw/standard_deviation_model_p_gaw

    #print normal_model_gaw
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
    #plt.plot(since2006_obs,var4, color='black', label='Mace Head Obs.')
    #plt.plot(since2006, model_cut, color='green', label='GEOS v9.01.03 Nested Europe ')
    #plt.grid(True)
    #leg=plt.legend(loc=first_label_pos)
    #leg.get_frame().set_alpha(0.4)
    #plt.xlabel('Time (Days since 2006)')
    #print units
    #plt.ylabel('%s (%s)' % (species_type,units))
    #plt.title('%s V Time' % (actual_species_name))

#Define sampling frequency
    samp_freq = 24

#Lomb-scargle plot
    ax= fig.add_subplot(1, 1, 1)

#Plot axis period lines and labels
    annotate_line_y=np.arange(1e-10,1e4,1)
    freq_year = [345]*len(annotate_line_y)
    plt.plot(freq_year, annotate_line_y,'r--',alpha=0.4)
    plt.text(345, 1e-10, '1 Year', fontweight='bold')

#Obs lomb
    fa, fb, nout, jmax, prob = lomb.fasper(since2006_obs, normal_var2, ofac, samp_freq)
#Divide output by sampling frequency
    fb = fb/samp_freq
    fb = np.log(fb)
    
    reversed_fb = fb[::-1]
    obs_periods = 1./fb
    reversed_obs_periods = obs_periods[::-1]
    obs_periods, obs_smoothed = ema(reversed_fb,reversed_obs_periods,150,20)
    #obs_smoothed=savitzky_golay(fb, window_size=301, order=1)
    obs_smoothed = np.exp(obs_smoothed) 

#Obs lomb
    fa_gaw, fb_gaw, nout, jmax, prob = lomb.fasper(since2006_obs_gaw, normal_var2_gaw, ofac, samp_freq)
#Divide output by sampling frequency
    fb_gaw = fb_gaw/samp_freq
    fb_gaw = np.log(fb_gaw)
   
    reversed_fb_gaw = fb_gaw[::-1]
    obs_periods_gaw = 1./fb_gaw
    reversed_obs_periods_gaw = obs_periods_gaw[::-1]
    obs_periods_gaw, obs_smoothed_gaw = ema(reversed_fb_gaw,reversed_obs_periods_gaw,150,20)
    #obs_smoothed_gaw=savitzky_golay(fb_gaw, window_size=301, order=1)
    obs_smoothed_gaw = np.exp(obs_smoothed_gaw)
    print 'obs_after_smooth', obs_smoothed_gaw
#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
    #nyquist_freq_lomb_obs = frequencies[-1]
    #Si_lomb_obs = np.mean(fb)*nyquist_freq_lomb_obs
    #print nyquist_freq_lomb_obs, Si_lomb_obs, Si_lomb_obs*2 

#plot up
    #plt.loglog(1./fa, fb,'kx',markersize=2, label='Mace Head Obs. ')

#Model lomb
    #print type(normal_model)
    fx, fy, nout, jmax, prob2 = lomb.fasper(since2006,normal_model, ofac, samp_freq)
#Divide output by sampling frequency
    fy = fy/samp_freq
    fy = np.log(fy)
   
    reversed_fy = fy[::-1]
    model_periods = 1./fx
    reversed_model_periods = model_periods[::-1]
    model_periods, model_smoothed = ema(reversed_fy,reversed_model_periods,150,20) 
    #model_smoothed=savitzky_golay(fy, window_size=301, order=1)
    model_smoothed = np.exp(model_smoothed)
 
#gaw basline lomb
    #print type(normal_model_gaw)
    #print normal_model_gaw
    fx_gaw, fy_gaw, nout, jmax, prob2 = lomb.fasper(since2006_gaw,normal_model_gaw, ofac, samp_freq)
#Divide output by sampling frequency
    fy_gaw = fy_gaw/samp_freq
    fy_gaw = np.log(fy_gaw)
   
    reversed_fy_gaw = fy_gaw[::-1]
    model_periods_gaw = 1./fx_gaw
    reversed_model_periods_gaw = model_periods_gaw[::-1]
    model_periods_gaw, model_smoothed_gaw = ema(reversed_fy_gaw,reversed_model_periods_gaw,150,20) 
    #model_smoothed_gaw=savitzky_golay(fy_gaw, window_size=301, order=1)
    model_smoothed_gaw = np.exp(model_smoothed_gaw)
    print 'model_after_smooth', model_smoothed_gaw
     #print model_smoothed_gaw

#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
    #nyquist_freq_lomb_model = frequencies[-1]
    #Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
    #print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 

#plot up
    #print obs_smoothed
    #print model_smoothed 
    #Which dataset is shorter
    # obs longer than model
    if len(obs_smoothed) > len(model_smoothed):
        print 'yes'
        obs_smoothed = obs_smoothed[:len(model_smoothed)]
        period_array = model_smoothed
	#period_array = model_periods[:len(model_smoothed)]
    #model longer than obs
    elif len(model_smoothed) >= len(obs_smoothed):
        print 'yes'
        model_smoothed = model_smoothed[:len(obs_smoothed)]
	period_array = obs_smoothed
        #period_array = obs_periods[:len(obs_smoothed)]

    if len(obs_smoothed_gaw) > len(model_smoothed_gaw):
        print 'yes'
        obs_smoothed_gaw = obs_smoothed_gaw[:len(model_smoothed_gaw)]
	period_array = model_periods_gaw
        #period_array_gaw = model_periods_gaw[:len(model_smoothed_gaw)]
    #model longer than obs
    elif len(model_smoothed_gaw) >= len(obs_smoothed_gaw):
        print 'gaw'
        model_smoothed_gaw = model_smoothed_gaw[:len(obs_smoothed_gaw)]
        #period_array_gaw = obs_periods_gaw[:len(obs_smoothed_gaw)]
        period_array = obs_periods_gaw

    

    print 'model_smoothed', model_smoothed_gaw
    print 'obs_smoothed', obs_smoothed_gaw
#calculate % of observations

        #covariance_array = np.hstack((fb,fy))

    compare_powers = model_smoothed/obs_smoothed
    compare_powers =  compare_powers *100

    compare_powers_gaw = model_smoothed_gaw/obs_smoothed_gaw
    compare_powers_gaw =  compare_powers_gaw *100

    #print compare_powers_gaw
   
    ax.set_xscale('log', basex=10)
    ax.set_yscale('log', basey=10)

    plt.plot(period_array, compare_powers , color='green', marker='x', alpha = 0.75, markersize=2, label = 'O3 % Diff Spatial correction.')
    plt.plot(period_array_gaw, compare_powers_gaw , color='black', marker='x', alpha = 0.3, markersize=2, label = 'O3 % Diff Baseline.')
    plt.grid(True)
    ax.xaxis.set_major_formatter(FormatStrFormatter('%.i'))
    ax.yaxis.set_major_formatter(FormatStrFormatter('%.i'))
    leg=plt.legend(loc=4, prop={'size':21})
    leg.get_frame().set_alpha(0.4)
    plt.xlim(0.05,1e1)
    plt.ylim(0.1,1000)
    plt.xlabel('Period (Days)', fontsize=21)
    plt.ylabel('Percent of Obs. PSD (%)', fontsize=21)
    plt.title('% PSD of Model compared to Obs.', fontsize=21)

    #percent1 = period_percent_diff(np.min(obs_periods),1,fb,fy,obs_periods,model_periods)
    #percent2 = period_percent_diff(1,2,fb,fy,obs_periods,model_periods)
    #percent3 = period_percent_diff(2,7,fb,fy,obs_periods,model_periods)
    
    #plt.grid(True)
    #leg=plt.legend(loc=7)
    #leg.get_frame().set_alpha(0.4)
    #plt.text(1e-2, 3000,'Period: 2 hours to 1 day, a %% Diff. of: %.2f%%'  %(percent1),  fontweight='bold')
    #plt.text(1e-2, 500,'Period: 1 day to 2 days, a %% Diff. of: %.2f%%'  %(percent2),  fontweight='bold')
    #plt.text(1e-2, 90,'Period: 2 days to 7 days, a %% Diff. of: %.2f%%'  %(percent3),  fontweight='bold')
    #plt.ylim(1e-10,1e4)
    #plt.xlabel('Period (Days)')
    #plt.ylabel(r'PSD $(ppb^{2}/days^{-1})$')
    #plt.title('Lomb-Scargle %s Power V Period' % actual_species_name)

#plt.savefig('O3_capeverde_comparison_plots.ps', dpi = 200)
    

    plt.show()



plot(species)
