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
def konnoOhmachiSmoothingWindow(frequencies, center_frequency, bandwidth,
                                normalize=True):
    """
    Returns the Konno & Ohmachi Smoothing window for every frequency in
    frequencies.

    Returns the smoothing window around the center frequency with one value per
    input frequency defined as follows (see [Konno1998]_):

    [sin(b * log_10(f/f_c)) / (b * log_10(f/f_c)]^4
        b   = bandwidth
        f   = frequency
        f_c = center frequency

    The bandwidth of the smoothing function is constant on a logarithmic scale.
    A small value will lead to a strong smoothing, while a large value of will
    lead to a low smoothing of the Fourier spectra.
    The default (and generally used) value for the bandwidth is 40. (From the
    Geopsy documentation - www.geopsy.org)

    All parameters need to be positive. This is not checked due to performance
    reasons and therefore any negative parameters might have unexpected
    results.

    This function might raise some numpy warnings due to divisions by zero and
    logarithms of zero. This is intentional and faster than prefiltering the
    special cases. You can disable numpy warnings (they usually do not show up
    anyways) with:

    temp = np.geterr()
    np.seterr(all='ignore')
    ...code that raises numpy warning due to division by zero...
    np.seterr(**temp)

    :param frequencies: numpy.ndarray (float32 or float64)
        All frequencies for which the smoothing window will be returned.
    :param center_frequency: float >= 0.0
        The frequency around which the smoothing is performed.
    :param bandwidth: float > 0.0
        Determines the width of the smoothing peak. Lower values result in a
        broader peak. Defaults to 40.
    :param normalize: boolean, optional
        The Konno-Ohmachi smoothing window is normalized on a logarithmic
        scale. Set this parameter to True to normalize it on a normal scale.
        Default to False.
    """
    if frequencies.dtype != np.float32 and frequencies.dtype != np.float64:
        msg = 'frequencies needs to have a dtype of float32/64.'
        raise ValueError(msg)
    # If the center_frequency is 0 return an array with zero everywhere except
    # at zero.
    if center_frequency == 0:
        smoothing_window = np.zeros(len(frequencies), dtype=frequencies.dtype)
        smoothing_window[frequencies == 0.0] = 1.0
        return smoothing_window
    # Calculate the bandwidth*log10(f/f_c)
    smoothing_window = bandwidth * np.log10(frequencies / center_frequency)
    # Just the Konno-Ohmachi formulae.
    smoothing_window[...] = (np.sin(smoothing_window) / smoothing_window) ** 4
    # Check if the center frequency is exactly part of the provided
    # frequencies. This will result in a division by 0. The limit of f->f_c is
    # one.
    smoothing_window[frequencies == center_frequency] = 1.0
    # Also a frequency of zero will result in a logarithm of -inf. The limit of
    # f->0 with f_c!=0 is zero.
    smoothing_window[frequencies == 0.0] = 0.0
    # Normalize to one if wished.
    if normalize:
        smoothing_window /= smoothing_window.sum()
    return smoothing_window


def calculateSmoothingMatrix(frequencies, bandwidth, normalize=True):
    """
    Calculates a len(frequencies) x len(frequencies) matrix with the Konno &
    Ohmachi window for each frequency as the center frequency.

    Any spectrum with the same frequency bins as this matrix can later be
    smoothed by a simple matrix multiplication with this matrix:
        smoothed_spectrum = np.dot(spectrum, smoothing_matrix)

    This also works for many spectra stored in one large matrix and is even
    more efficient.

    This makes it very efficient for smoothing the same spectra again and again
    but it comes with a high memory consumption for larger frequency arrays!

    :param frequencies: numpy.ndarray (float32 or float64)
        The input frequencies.
    :param bandwidth: float > 0.0
        Determines the width of the smoothing peak. Lower values result in a
        broader peak. Defaults to 40.
    :param normalize: boolean, optional
        The Konno-Ohmachi smoothing window is normalized on a logarithmic
        scale. Set this parameter to True to normalize it on a normal scale.
        Default to False.
    """
    # Create matrix to be filled with smoothing entries.
    sm_matrix = np.empty((len(frequencies), len(frequencies)),
                         frequencies.dtype)
    for _i, freq in enumerate(frequencies):
        sm_matrix[_i, :] = konnoOhmachiSmoothingWindow(frequencies, freq,
                                              bandwidth, normalize=normalize)
    return sm_matrix

def konnoOhmachiSmoothing(spectra, frequencies, bandwidth=40, count=1,
                  enforce_no_matrix=False, max_memory_usage=512,
                  normalize=True):
    """
    Smoothes a matrix containing one spectra per row with the Konno-Ohmachi
    smoothing window.

    All spectra need to have frequency bins corresponding to the same
    frequencies.

    This method first will estimate the memory usage and then either use a fast
    and memory intensive method or a slow one with a better memory usage.

    :param spectra: numpy.ndarray (float32 or float64)
        One or more spectra per row. If more than one the first spectrum has to
        be accessible via spectra[0], the next via spectra[1], ...
    :param frequencies: numpy.ndarray (float32 or float64)
        Contains the frequencies for the spectra.
    :param bandwidth: float > 0.0
        Determines the width of the smoothing peak. Lower values result in a
        broader peak. Defaults to 40.
    :param count: integer, optional
        How often the apply the filter. For very noisy spectra it is useful to
        apply is more than once. Defaults to 1.
    :param enforce_no_matrix: boolean, optional
        An efficient but memory intensive matrix-multiplication algorithm is
        used in case more than one spectra is to be smoothed or one spectrum is
        to be smoothed more than once if enough memory is available. This flag
        disables the matrix algorithm altogether. Defaults to False
    :param max_memory_usage: integer, optional
        Set the maximum amount of extra memory in MB for this method. Decides
        whether or not the matrix multiplication method is used. Defaults to
        512 MB.
    :param normalize: boolean, optional
        The Konno-Ohmachi smoothing window is normalized on a logarithmic
        scale. Set this parameter to True to normalize it on a normal scale.
        Default to False.
    """
    if (frequencies.dtype != np.float32 and frequencies.dtype != np.float64) \
       or (spectra.dtype != np.float32 and spectra.dtype != np.float64):
        msg = 'frequencies and spectra need to have a dtype of float32/64.'
        raise ValueError(msg)
    # Spectra and frequencies should have the same dtype.
    if frequencies.dtype != spectra.dtype:
        frequencies = np.require(frequencies, np.float64)
        spectra = np.require(spectra, np.float64)
        msg = 'frequencies and spectra should have the same dtype. It ' + \
              'will be changed to np.float64 for both.'
        warnings.warn(msg)
    # Check the dtype to get the correct size.
    if frequencies.dtype == np.float32:
        size = 4.0
    elif frequencies.dtype == np.float64:
        size = 8.0
    # Calculate the approximate usage needs for the smoothing matrix algorithm.
    #length = len(frequencies)
    #approx_mem_usage = 100 #(length * length + 2 * len(spectra) + length) * \
            #size / 1048576.0
    # If smaller than the allowed maximum memory consumption build a smoothing
    # matrix and apply to each spectrum. Also only use when more then one
    # spectrum is to be smoothed.
    #if enforce_no_matrix is False and (len(spectra.shape) > 1 or count > 1) \
    #   and approx_mem_usage < max_memory_usage:
        # Disable numpy warnings due to possible divisions by zero/logarithms
        # of zero.
    #if 1 ==1:
     #   print 'smaller'
     #   temp = np.geterr()
     #   np.seterr(all='ignore')
     #   smoothing_matrix = calculateSmoothingMatrix(frequencies, bandwidth,
     #                                        normalize=normalize)
     #   np.seterr(**temp)
     #   new_spec = np.dot(spectra, smoothing_matrix)
       # Eventually apply more than once.
     #   for _i in xrange(count - 1):
     #       new_spec = np.dot(new_spec, smoothing_matrix)
     #   return new_spec
    # Otherwise just calculate the smoothing window every time and apply it.
    if 1==1:
        new_spec = np.empty(spectra.shape, spectra.dtype)
        print 'more'
        # Separate case for just one spectrum.
        if len(new_spec.shape) == 1:
            # Disable numpy warnings due to possible divisions by
            # zero/logarithms of zero.
            temp = np.geterr()
            np.seterr(all='ignore')
            for _i in xrange(len(frequencies)):
                window = konnoOhmachiSmoothingWindow(frequencies,
                        frequencies[_i], bandwidth, normalize=normalize)
                new_spec[_i] = (window * spectra).sum()
            np.seterr(**temp)
        # Reuse smoothing window if more than one spectrum.
        else:
            # Disable numpy warnings due to possible divisions by
            # zero/logarithms of zero.
            temp = np.geterr()
            np.seterr(all='ignore')
            for _i in xrange(len(frequencies)):
                window = konnoOhmachiSmoothingWindow(frequencies,
                        frequencies[_i], bandwidth, normalize=normalize)
                for _j, spec in enumerate(spectra):
                    new_spec[_j, _i] = (window * spec).sum()
            np.seterr(**temp)
        # Eventually apply more than once.
        while count > 1:
            new_spec = konnoOhmachiSmoothing(new_spec, frequencies, bandwidth,
                                enforce_no_matrix=True, normalize=normalize)
            count -= 1
        return new_spec





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
def plot():

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


#now read in the observations

    myfile=nappy.openNAFile('York_merge_Cape_verde_1hr_R1.na')
    myfile.readData()

#ppy.openNAFile('York_merge_Cape_verde_1hr_R1.na')
    counter = 0
    fig =plt.figure(figsize=(20,12)) 
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
            model_cut = model[:,k]*unit_cut
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


#Define sampling frequency
        samp_freq = 24

#Lomb-scargle plot

#Plot axis period lines and labels
        annotate_line_y=np.arange(1e-10,1e4,1)
        horiz_line_100 =np.arange(0,2000,1)
        freq_year = [345]*len(annotate_line_y)
        array_100 = [100]*len(horiz_line_100)
        plt.plot(freq_year, annotate_line_y,'r--',alpha=0.4)
        plt.text(345, 5, '1 Year', fontweight='bold')
        plt.plot(horiz_line_100, array_100,'r--',alpha=0.4)
        plt.text(1024, 80, '100%', fontweight='bold')

#Obs lomb
        fa, fb, nout, jmax, prob = lomb.fasper(time2, normal_var2, ofac, samp_freq)
        obs_sig = fa, fb, nout, ofac
#Divide output by sampling frequency
        fb = fb/samp_freq

        len_fb = len(fb)

        zeropad = np.zeros(10000)
        fb = np.concatenate((fb,zeropad))
        padded_obs_period = np.concatenate((fa,zeropad))
        
        obs_smoothed = konnoOhmachiSmoothing(fb, padded_obs_period, bandwidth=40, count=1,
                  enforce_no_matrix=True, max_memory_usage=512,
                  normalize=False)

        obs_smoothed = obs_smoothed[:len_fb]
#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
    #nyquist_freq_lomb_obs = frequencies[-1]
    #Si_lomb_obs = np.mean(fb)*nyquist_freq_lomb_obs
    #print nyquist_freq_lomb_obs, Si_lomb_obs, Si_lomb_obs*2 

#plot up
    #plt.loglog(1./fa, fb,'kx',markersize=2, label='Cape Verde Obs. ')

#Model lomb
        fx, fy, nout, jmax, prob2 = lomb.fasper(since2006,normal_model, ofac, samp_freq)
        model_sig = fx, fy, nout, ofac
#Divide output by sampling frequency
        fy = fy/samp_freq
        
        len_fy = len(fy)    
   
        fy = np.concatenate((fy,zeropad))
        padded_model_period = np.concatenate((fx,zeropad))
     
        model_smoothed = konnoOhmachiSmoothing(fy, padded_model_period, bandwidth=40, count=1,
                  enforce_no_matrix=True, max_memory_usage=512,
                  normalize=False)


        model_smoothed = model_smoothed[:len_fy]
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
            freq_array = fa
            period_array = obs_periods

#calculate % of observations
    
        #covariance_array = np.hstack((fb,fy))
        
        compare_powers = model_smoothed/obs_smoothed 
        compare_powers =  compare_powers *100  
        #compare_powers =  np.cov(covariance_array, y=None, rowvar=0, bias=0, ddof=None)
        #print compare_powers 
        #window_size = 120
        #compare_powers = np.log(compare_powers)
        #compare_powers = movingaverage(compare_powers,window_size)
        #compare_powers = smooth(compare_powers,window_size)
        #compare_powers = 10**compare_powers


        #cut_powers = [y for x,y in zip(obs_periods,compare_powers) if x < end]
        #cut_periods = [x for x,y in zip(obs_periods,compare_powers) if x < end]
        #rest_powers = [y for x,y in zip(obs_periods,compare_powers) if  x >= end]
        #rest_cut_periods = [x for x,y in zip(obs_periods,compare_powers) if x >= end]
        #smooth(compare_powers,window_size,obs_periods,species,counter)
        #for i in rest_powers:
        #    compare_powers=np.append(compare_powers, i)
        #print compare_powers
        #for i in rest_cut_periods:
        #    cut_periods=np.append(cut_periods,i)
        #compare_powers = compare_powers.flatten() 
        #compare_periods = combined_periods.flatten() 
        #cut_periods = cut_periods.flatten() 
        #print rest_cut_periods[-1]
        #print combined_periods[-1]
       

        #smoothed = konnoOhmachiSmoothing(compare_powers, freq_array, bandwidth=40, count=1,
        #          enforce_no_matrix=True, max_memory_usage=512,
        #          normalize=False) 
        
        # xticks = [0.08,0.10,0.12,0.15,0.25,0.50,1,10,100,1000,5000]
        #ax.set_xticks(xticks)


        #ax.set_xticklabels(xticks)
        #locator = LinearLocator 
        #ax.xaxis.set_major_locator(locator)

        #standard_deviation_analysis = np.std(compare_powers)
        #mean_analysis = np.mean(compare_powers)
        #normal_analysis = compare_powers-mean_analysis
        #normal_analysis = normal_analysis/standard_deviation_analysis

        ax.set_xscale('log', basex=10)
        ax.set_yscale('log', basey=10)
       
          
        plt.plot(period_array, compare_powers , color=colour_list[counter], marker='x', alpha = 0.75, markersize=2, label = species)
        #ax.plot(rest_cut_periods, rest_powers , color=colour_list[counter], marker='x', alpha = 0.75, markersize=2, label = species)  
  #percent1 = period_percent_diff(np.min(obs_periods),1,fb,fy,obs_periods,model_periods)
    #percent2 = period_percent_diff(1,2,fb,fy,obs_periods,model_periods)
    #percent3 = period_percent_diff(2,7,fb,fy,obs_periods,model_periods)
    
        plt.grid(True)
        ax.xaxis.set_major_formatter(FormatStrFormatter('%.i'))
        ax.yaxis.set_major_formatter(FormatStrFormatter('%.i'))
        leg=plt.legend(loc=4)
        leg.get_frame().set_alpha(0.4)
    #plt.text(1e-2, 3000,'Period: 2 hours to 1 day, a %% Diff. of: %.2f%%'  %(percent1),  fontweight='bold')
    #plt.text(1e-2, 500,'Period: 1 day to 2 days, a %% Diff. of: %.2f%%'  %(percent2),  fontweight='bold')
    #plt.text(1e-2, 90,'Period: 2 days to 7 days, a %% Diff. of: %.2f%%'  %(percent3),  fontweight='bold')
        plt.ylim(0.05,10000)
        plt.xlabel('Period (Days)')
        plt.ylabel('Percent of Obs. PSD (%)')
        plt.title('% PSD of Model compared to Obs.')
        counter+=1
#plt.savefig('O3_capeverde_comparison_plots.ps', dpi = 200)
    
     
    plt.show()


plot()
