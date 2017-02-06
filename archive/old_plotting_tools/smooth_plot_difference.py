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


    #make index for high frequency measure of obs. and model PSD. Say is between min period(2 hours) and 1 day 
def period_percent_diff(start, end, fb, fy, obs_periods, model_periods):
        
    cut_obs_periods = [x for x in obs_periods if start <= x < end]
    cut_model_periods = [x for x in model_periods if start <= x < end]

    #global fb
    #global fy
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



def readfile(filename, location):

    read = np.load(filename)
    locs = np.where(read[:,0] == location)
    big = read[locs]
    big = np.float64(big)

    return big


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
    for files in filename:
        print files
        reader=csv.reader(open(files,'rb'), delimiter=',', skipinitialspace=True)
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
            print i
        vals.append(i)
    vals = np.array(vals)
    vals = vals.astype(np.float)
    #vals = vals.astype(int)
    return date, time, vals

#reads in species name and processes unique variables
def plot(species,location):
   
    #Set obs data for each location

    if location == 'Arrival_Heights':
        obsfile = 'arrival_heights_o3_hourly/o3*'
        loc_label = 'Arrival Heights'
        gaw_code = 010  
 
    elif location == 'Barrow':
        obsfile = 'barrow_o3_hourly/o3*'
        loc_label = 'Barrow'
        gaw_code = 015  
   
    elif location == 'Lauder':
        obsfile = 'lauder_o3_hourly/o3*'
        loc_label = 'Lauder'
        gaw_code = 106 
 
    elif location == 'Mace_Head':
        obsfile = 'O3_mace_head_ppbV.txt'
        loc_label = 'Mace Head'    
        gaw_code = 112    

    elif location == 'Mauna_Loa':
        obsfile = 'mauna_loa_o3_hourly/o3*'
        loc_label = 'Mauna Loa'
        gaw_code = 116 

    elif location == 'Niwot_Ridge':
        obsfile = 'niwot_ridge_o3_hourly/o3*'
        loc_label = 'Niwot Ridge'
        gaw_code = 132 

    elif location == 'Ragged_Point':
        obsfile = 'ragged_point_o3_hourly/o3*'
        loc_label = 'Ragged Point'
        gaw_code = 152 
 
    elif location == 'South_Pole':
        obsfile = 'south_pole_o3_hourly/o3*'
        loc_label = 'South Pole'
        gaw_code = 173 

    elif location == 'Trinidad_Head':
        obsfile = 'trinidad_head_o3_hourly/o3*'
        loc_label = 'Trinidad Head'
        gaw_code = 189 

    elif location == 'Tudor_Hill':
        obsfile = 'tudor_hill_o3_hourly/o3*'
        loc_label = 'Tudor Hill'
        gaw_code = 191 

    elif location == 'Tutuila':
        obsfile = 'tutuila_o3_hourly/o3*'
        loc_label = 'Tutuila'
        gaw_code = 192 

    #Set model_cut switch to default 0, if want to do more complicated cuts from model field, specify model_cut_switch == 1 in species definitions
    model_cut_switch = 0
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
        first_label_pos = 1
        obs_data_name = 'Air Temperature (degC) Campbell_(Mean)'
        unit_cut= 1
        species_type = 'Temp.'
        actual_species_name = 'Surface Temperature'

    elif species == 'GMAO_PSFC': #Surface Pressure
        units = 'hPa'
        first_label_pos = 1
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
        first_label_pos = 1
        obs_data_name = 'Wind Speed (m/s) Campbell_(Mean)'
        unit_cut= 1
        species_type = 'Wind Speed'
        model_cut_switch = 1
        actual_species_name = 'Surface Windspeed'

    elif species == 'GMAO_RADSW': #Sensible heat flux form surface       
        units = r'$Wm^{-2}$'
        first_label_pos = 1
        obs_data_name = 'Solar Radiation (Wm-2) Campbell_(Mean)'
        unit_cut= 1
        species_type = 'Solar Radiation'
        actual_species_name = 'Surface Solar Radiation'

    elif species == 'GMAO_ABSH': #Absolute Humidity       
        units = 'molec/cm-3'
        first_label_pos = 1
        obs_data_name = ''
        unit_cut= 1
        species_type = 'Absolute Humidity'
        actual_species_name = 'Absolute Humidity'

    elif species == 'GMAO_RHUM': #Relative Humidity       
        units = '%'
        first_label_pos = 1
        obs_data_name = 'Relative Humidity (%) Campbell_(Mean)'
        unit_cut= 1
        species_type = 'Relative Humidity'
        actual_species_name = 'Relative Humidity'
    
   
#do I need to read everything in 
    try: 
        names
    except NameError:
# Readin the model output

        model = readfile("gaw_logs.npy", gaw_code)
       
        print model.shape
        print model
 
# Processes the date 
        year=(model[:,1]//10000)
        month=((model[:,1]-year*10000)//100)
        day=(model[:,1]-year*10000-month*100)
    
        hour=model[:,2]//100
        min=(model[:,2]-hour*100)
    
        doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(2006,1,1,0,0,0) \
              for i in range(len(year))]
    
        since2006=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]


#now read in the observations
    if location == 'Mace_Head': 
        date, time, vals =  NOAA_data_reader_mace_head(glob.glob(obsfile))
    else:
        date, time, vals =  NOAA_data_reader(glob.glob(obsfile))
    valid = vals >= 0
    vals = vals[valid]
    date = date[valid]
    time = time[valid]
    print vals
# Process NOAA obs time
    year=(date//10000)
    month=((date-year*10000)//100)
    day=(date-year*10000-month*100)

    hour=time//100
    min=(time-hour*100)

    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(2006,1,1,0,0,0) \
              for i in range(len(year))]

    since2006_2=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]

#Pre normalise obs data for lomb analysis
    standard_deviation_obs_p = np.std(vals)
    mean_obs_p = np.mean(vals)
    normal_var2 = vals-mean_obs_p
    normal_var2 = normal_var2/standard_deviation_obs_p
 
#Calculate variance of pre-processed obs data- should be 1 if normal
    #standard_dev_obs = np.std(normal_var_2, dtype=np.float64)
    #variance_obs = standard_dev_obs**2
    #print 'Variance - pre-processed obs data= ', variance_obs

#Define sampling intervals
    samp_spacing = 1./24.

#Convert model time array into numpy array
    since2006=np.array(since2006)
    since2006_2=np.array(since2006_2)
#Need to normalise model data also
    model_cut = model[:,3]*unit_cut  	
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
    ax = plt.subplot(111)
#Plot up standard conc. v time plot
    #ax1= fig.add_subplot(2, 1, 1)
    #fig.subplots_adjust(hspace=0.3)
    #plt.plot(since2006_2,vals, color='black', label= '%s Obs.' % loc_label)
    #plt.plot(since2006, model_cut, color='green', label='GEOS v9.01.03 4x5 ')
    #plt.grid(True)
    #leg=plt.legend(loc=1)
    #leg.get_frame().set_alpha(0.4)
    #plt.xlabel('Time (Days since 2006)')
    #print units
    #plt.ylabel('%s (%s)' % (species_type,units))
    #plt.title('%s V Time' % (actual_species_name))

#Define sampling frequency
    samp_freq = 24

#Lomb-scargle plot
    #ax3= fig.add_subplot(2, 1, 2)

#Plot axis period lines and labels
    annotate_line_y=np.arange(1e-10,1e4,1)
    horiz_line_100 =np.arange(0,2000,1)
    freq_year = [345]*len(annotate_line_y)
    array_100 = [100]*len(horiz_line_100)
    plt.plot(freq_year, annotate_line_y,'r--',alpha=0.4)
    plt.text(345, 10, '1 Year', fontweight='bold')
    #plt.plot(horiz_line_100, array_100,'r--',alpha=0.4)
    #plt.text(0.05, 80, '100%', fontweight='bold')

#Obs lomb
    fa, fb, nout, jmax, prob = lomb.fasper(since2006_2, normal_var2, ofac, samp_freq)
#Divide output by sampling frequency
    fb = fb/samp_freq
    fb = np.log(fb)
    obs_smoothed=savitzky_golay(fb, window_size=301, order=1)
    obs_smoothed = np.exp(obs_smoothed) 

    #nyquist_freq_lomb_obs = frequencies[-1]
    #Si_lomb_obs = np.mean(fb)*nyquist_freq_lomb_obs
    #print nyquist_freq_lomb_obs, Si_lomb_obs, Si_lomb_obs*2 

#plot up
    #plt.loglog(1./fa, fb,'kx',markersize=2, label='%s Obs. ' % loc_label)
#Model lomb
    #print normal_model
    fx, fy, nout, jmax, prob2 = lomb.fasper(since2006,normal_model, ofac, samp_freq)
#Divide output by sampling frequency
    fy = fy/samp_freq
    fy = np.log(fy)
    model_smoothed=savitzky_golay(fy, window_size=301, order=1)
    model_smoothed = np.exp(model_smoothed)
    #print model_smoothed
#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
    #nyquist_freq_lomb_model = frequencies[-1]
    #Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
    #print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 

#plot up
    #plt.loglog(1./fx, fy, 'gx', alpha = 0.75,markersize=2, label='GEOS v9.01.03 4x5 ')
    #plt.loglog(1./fa, obs_smoothed, color = 'orangered', marker='x',linestyle='None',  alpha = 0.75,markersize=2, label='Smoothed Mace Head Obs. ')    
    #plt.loglog(1./fx, model_smoothed, color = 'blue', marker='x', linestyle='None',  alpha = 0.75,markersize=2, label='Smoothed GEOS v9.01.03 4x5 ')    
    
    obs_periods = 1./fa
    model_periods = 1./fx

     
    #Which dataset is shorter
    # obs longer than model
    if len(obs_smoothed) > len(model_smoothed):
        obs_smoothed = obs_smoothed[:len(model_smoothed)]
        freq_array = fx
        period_array = model_periods
    #model longer than obs
    if len(model_smoothed) > len(obs_smoothed):
        model_smoothed = model_smoothed[:len(obs_smoothed)]
        freq_array = fa
        period_array = obs_periods

#calculate % of observations

        #covariance_array = np.hstack((fb,fy))
   
   # print model_smoothed 
    compare_powers = model_smoothed/obs_smoothed
    compare_powers =  compare_powers *100
   
    ax.set_xscale('log', basex=10)
    ax.set_yscale('log', basey=10)


    plt.plot(period_array, compare_powers , color='black', marker='x', alpha = 0.75, markersize=2, label = 'O3 % Diff.')  
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
 

#plt.savefig('O3_capeverde_comparison_plots.ps', dpi = 200)
    plt.show()


