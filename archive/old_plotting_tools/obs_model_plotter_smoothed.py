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

species_list = ['O3']#,'CO','C2H6','C3H8','TEMP','SURFACE_PRES','WINDSPEED','REL_HUMIDITY']
colour_list=['LimeGreen','Red','Blue','yellow','Magenta','Cyan','DarkOrange','black','SteelBlue']
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
#   2     r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
#   3     The Savitzky-Golay filter removes high frequency noise from data.
#   4     It has the advantage of preserving the original shape and
#   5     features of the signal better than other types of filtering
#   6     approaches, such as moving averages techniques.
#   7     Parameters
#   8     ----------
#   9     y : array_like, shape (N,)
#  10         the values of the time history of the signal.
#  11     window_size : int
#  12         the length of the window. Must be an odd integer number.
#  13     order : int
#  14         the order of the polynomial used in the filtering.
#  15         Must be less then `window_size` - 1.
#  16     deriv: int
#  17         the order of the derivative to compute (default = 0 means only smoothing)
#  18     Returns
#  19     -------
#  20     ys : ndarray, shape (N)
#  21         the smoothed signal (or it's n-th derivative).
#  22     Notes
#  23     -----
#  24     The Savitzky-Golay is a type of low-pass filter, particularly
#  25     suited for smoothing noisy data. The main idea behind this
#  26     approach is to make for each point a least-square fit with a
#  27     polynomial of high order over a odd-sized window centered at
#  28     the point.
#  29     Examples
#  30     --------
#  31     t = np.linspace(-4, 4, 500)
#  32     y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
#  33     ysg = savitzky_golay(y, window_size=31, order=4)
#  34     import matplotlib.pyplot as plt
#  35     plt.plot(t, y, label='Noisy signal')
#  36     plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
#  37     plt.plot(t, ysg, 'r', label='Filtered signal')
#  38     plt.legend()
#  39     plt.show()
#  40     References
#  41     ----------
#  42     .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
#  43        Data by Simplified Least Squares Procedures. Analytical
#  44        Chemistry, 1964, 36 (8), pp 1627-1639.
#  45     .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
#  46        W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
#  47        Cambridge University Press ISBN-13: 978052188068"""
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
    names = read[0,2:]
    locs = np.where(read[:,1] == location)
    big = read[locs]
    valid_data = big[:,2:]
    names = names.tolist()
    valid_data = np.float64(valid_data)

    return valid_data, names

def smoother(x,window_len):
    s=np.r_[2*x[0]-x[window_len-1::-1],x,2*x[-1]-x[-1:-window_len:-1]]
    w=np.ones(window_len,'d')
    y=np.convolve(w/w.sum(),s,mode='same')
    calc = y[window_len:-window_len+1]
    return calc



def smooth(x,window_len,obs_periods,species,counter):
    #n_intervals =  len(x) // window_len  
    
    sorted_periods, sorted_powers = zip(*sorted(zip(obs_periods,x)))
    #sorted_periods, sorted_powers = (list(i) for i in zip(*sorted(zip(obs_periods, x))))
    sorted_periods= list(sorted_periods)
    sorted_powers= list(sorted_powers)
    zones=20

    remainder = np.mod(len(x), zones)
    remainder_powers = sorted_powers[-remainder:] 
    remainder_periods = sorted_periods[-remainder:]
   
    new_sorted_powers = sorted_powers[:-remainder]
    new_sorted_periods = sorted_periods[:-remainder]    

    n_intervals =  len(x) // zones 
    # Split into 10 chunks
    #reshaped_powers = np.reshape(new_sorted_powers,(-1,window_len))
    #reshaped_periods = np.reshape(new_sorted_periods,(-1,window_len))
    
    reshaped_powers = np.reshape(new_sorted_powers,(-1,n_intervals))
    reshaped_periods = np.reshape(new_sorted_periods,(-1,n_intervals))  
    window_len = 1000
    count = 0
    for i,j in zip(reshaped_powers,reshaped_periods):
        if window_len > 1:
            #start = 0 + n_intervals
            #stop = window_len+count
            #i = new_sorted_powers[start:stop]
            #j = new_sorted_periods[start:stop]
            count+=1
            s=np.r_[2*i[0]-i[window_len-1::-1],i,2*i[-1]-i[-1:-window_len:-1]]
            w=np.ones(window_len,'d')
            y=np.convolve(w/w.sum(),s,mode='same')
            calc = y[window_len:-window_len+1]
            plt.plot(j, calc , color=colour_list[counter], marker='x', alpha = 0.75, markersize=2, label = species)
            window_len= window_len-50
        else:
            plt.plot(j, calc , color=colour_list[counter], marker='x', alpha = 0.75, markersize=2, label = species)
    i = remainder_powers
    window_len = remainder
    print window_len
    s=np.r_[2*i[0]-i[window_len-1::-1],i,2*i[-1]-i[-1:-window_len:-1]]
    w=np.ones(window_len,'d')
    y=np.convolve(w/w.sum(),s,mode='same')
    calc = y[window_len:-window_len+1]
    plt.plot(remainder_periods,  calc , color=colour_list[counter], marker='x', alpha = 0.75, markersize=2, label = species)
#def movingaverage(interval, window_size):
#    n_intervals =  len(interval) // window_size  
#    remainder = np.mod(len(interval), window_size)
#    remainder_array = interval[-remainder:] 
#    new_interval = interval[:-remainder]

#    reshaped_interval = np.reshape(new_interval,(-1,n_intervals))
#    big_array=[]
#    for i in reshaped_interval:
#        if len(i) < window_size:
#            window_size = len(i)
        
        #big_window= np.ones(int(window_size))/float(window_size)
        #calc = np.convolve(i, big_window, 'same')
#        big_array.append(calc)
#        window_size= window_size-1
         
 
#    big_array = np.array(big_array)
#    big_array = big_array.flatten()
#    array = np.concatenate((big_array,remainder_array))
#    return array

#make index for high frequency measure of obs. and model PSD. Say is between min period(2 hours) and 1 day 
                                     
#reads in species name and processes unique variables
def plot():

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


#now read in the observations

    myfile=nappy.openNAFile('York_merge_Cape_verde_1hr_R1.na')
    myfile.readData()

#ppy.openNAFile('York_merge_Cape_verde_1hr_R1.na')
    counter = 0
    fig =plt.figure(figsize=(20,12)) 
    fig.patch.set_facecolor('white')
    ax = plt.subplot(111)
    
    for species in species_list:
    #Gives species exact model tags for convenience
        print species
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

        model_cut_switch = 0
        obs_switch = 0
        ofac = 1
        if species == 'O3':
            print 'yes'
            Units = 'ppbV'
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
        since2006=np.array(since2006)

#Need to normalise model data also
        if model_cut_switch == 0:
            k=names.index(species)
            print model[:,k]
            model_cut = model[:,k]*unit_cut
        if model_cut_switch == 1:
            model_cut = read_diff_species()
        

   #Add seasonal emission trend onto ethane.
        first_season = np.linspace(0,100,num=91,endpoint=True)
        second_season = first_season[::-1]
        third_season = np.linspace(0,-100, num=91, endpoint=True)
        fourth_season = third_season[::-1]
        fourth_season =np.append(fourth_season,0)


        n=0
        n_end=24
        step = 24
        season_index = 0
        new_model_cut=[]
        year_count = 0

        count=0
        while 1==1:
            if count <91:
                sliced = model_cut[n:n_end]
                season_value = first_season[season_index]
                sliced = [a+season_value for a in sliced]
                new_model_cut.append(sliced)
                n+=step
                n_end+=step
                #print 'season_1', count, season_index
                season_index+=1
                if season_index == 91:
                    season_index= 0

            elif count <182:
                sliced = model_cut[n:n_end]
                season_value = second_season[season_index]
                sliced = [a+season_value for a in sliced]
                new_model_cut.append(sliced)
                n+=step
                n_end+=step
                #print 'season_2', count, season_index
                season_index+=1
                if season_index == 91:
                    season_index= 0

            elif count < 273:
                sliced = model_cut[n:n_end]
                season_value = third_season[season_index]
                sliced = [a+season_value for a in sliced]
                new_model_cut.append(sliced)
                n+=step
                n_end+=step
                #print 'season_3', count, season_index
                season_index+=1
                if season_index == 91:
                    season_index= 0

            elif count < 365:
                sliced = model_cut[n:n_end]
                season_value = fourth_season[season_index]
                sliced = [a+season_value for a in sliced]
                new_model_cut.append(sliced)
                n+=step
                n_end+=step
                #print 'season_4', count, season_index
                season_index+=1

            else:
                count = 0
                year_count+=1
                season_index = 0
                continue

            if year_count == 6:
                break

            count+=1

        new_model_cut = reduce(lambda x,y: x+y,new_model_cut)

        standard_deviation_model_p = np.std(model_cut)
        mean_model_p = np.mean(model_cut)
        normal_model = model_cut-mean_model_p
        normal_model = normal_model/standard_deviation_model_p

        standard_deviation_model_p_corrected = np.std(new_model_cut)
        mean_model_p_corrected = np.mean(new_model_cut)
        normal_model_corrected = new_model_cut-mean_model_p_corrected
        normal_model_corrected = normal_model_corrected/standard_deviation_model_p_corrected

#Calculate variance of pre-processed model data- should be 1 if normal
    #standard_dev_model = np.std(normal_model, dtype=np.float64)
    #variance_model = standard_dev_model**2
    #print 'Variance - pre-processed model data= ', variance_model


#Define sampling frequency
        samp_freq = 24

#Lomb-scargle plot

#Plot axis period lines and labels
        #annotate_line_y=np.arange(1e-10,1e4,1)
        #horiz_line_100 =np.arange(0,2000,1)
        #freq_year = [345]*len(annotate_line_y)
        #array_100 = [100]*len(horiz_line_100)
        #plt.plot(freq_year, annotate_line_y,'r--',alpha=0.4)
        #plt.text(345, 5, '1 Year', fontweight='bold')
        #plt.plot(horiz_line_100, array_100,'r--',alpha=0.4)
        #plt.text(1024, 80, '100%', fontweight='bold')

#Obs lomb
        fa, fb, nout, jmax, prob = lomb.fasper(time2, normal_var2, ofac, samp_freq)
#Divide output by sampling frequency
        fb = fb/samp_freq
        fb = np.log(fb)
        obs_smoothed=savitzky_golay(fb, window_size=301, order=1)
        obs_smoothed = np.exp(obs_smoothed)
   
#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
    #nyquist_freq_lomb_obs = frequencies[-1]
    #Si_lomb_obs = np.mean(fb)*nyquist_freq_lomb_obs
    #print nyquist_freq_lomb_obs, Si_lomb_obs, Si_lomb_obs*2 

#plot up
    #plt.loglog(1./fa, fb,'kx',markersize=2, label='Cape Verde Obs. ')

#Model lomb
        fx, fy, nout, jmax, prob2 = lomb.fasper(since2006,normal_model, ofac, samp_freq)
#Divide output by sampling frequency
        fy = fy/samp_freq
        fy = np.log(fy) 
        model_smoothed=savitzky_golay(fy, window_size=301, order=1)   
        model_smoothed = np.exp(model_smoothed) 

#Model lomb
        fx, fy, nout, jmax, prob2 = lomb.fasper(since2006,normal_model_corrected, ofac, samp_freq)
#Divide output by sampling frequency
        fy_corrected = fy/samp_freq
        fy_corrected = np.log(fy)
        model_corrected_smoothed=savitzky_golay(fy_corrected, window_size=301, order=1)
        model_corrected_smoothed = np.exp(model_corrected_smoothed)
#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
    #nyquist_freq_lomb_model = frequencies[-1]
    #Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
    #print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 

#plot up
    #plt.loglog(1./fx, fy, 'gx', alpha = 0.75,markersize=2, label='GEOS v9.01.03 4x5 ')

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
            model_corrected_smoothed =  model_corrected_smoothed[:len(obs_smoothed)]
            freq_array = fa
            period_array = obs_periods

#calculate % of observations
    
        #covariance_array = np.hstack((fb,fy))
        
        compare_powers = model_smoothed/obs_smoothed 
        compare_powers =  compare_powers *100

        corrected_compare_powers = model_corrected_smoothed/obs_smoothed
        corrected_compare_powers =  corrected_compare_powers *100  

        ax.set_xscale('log', basex=10)
        ax.set_yscale('log', basey=10)
       
        #plt.plot(obs_periods,fb, color = 'k', marker='x', alpha = 0.75, markersize=2, label = 'Mace Head'  
       
        #plt.plot(period_array, corrected_compare_powers , color=colour_list[counter], marker='x', alpha = 0.75, markersize=2, label = species)
        plt.plot(period_array, compare_powers , color='black', marker='x', alpha = 0.75, markersize=2, label = species)
        #ax.plot(rest_cut_periods, rest_powers , color=colour_list[counter], marker='x', alpha = 0.75, markersize=2, label = species)  
  #percent1 = period_percent_diff(np.min(obs_periods),1,fb,fy,obs_periods,model_periods)
    #percent2 = period_percent_diff(1,2,fb,fy,obs_periods,model_periods)
    #percent3 = period_percent_diff(2,7,fb,fy,obs_periods,model_periods)
    
        plt.grid(True)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.i'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.i'))
        leg=plt.legend(loc=4, prop={'size':21})
        leg.get_frame().set_alpha(0.4)
    #plt.text(1e-2, 3000,'Period: 2 hours to 1 day, a %% Diff. of: %.2f%%'  %(percent1),  fontweight='bold')
    #plt.text(1e-2, 500,'Period: 1 day to 2 days, a %% Diff. of: %.2f%%'  %(percent2),  fontweight='bold')
    #plt.text(1e-2, 90,'Period: 2 days to 7 days, a %% Diff. of: %.2f%%'  %(percent3),  fontweight='bold')
        plt.xlim(0.05,1e1)
        plt.ylim(0.001,1e3)
        plt.xlabel('Period (Days)', fontsize=21)
        plt.ylabel('Percent of Obs. PSD (%)', fontsize=21)
        plt.title('% PSD of Model compared to Obs.',fontsize=21)
        counter+=1
#plt.savefig('O3_capeverde_comparison_plots.ps', dpi = 200)
    
     
    plt.show()


plot()
