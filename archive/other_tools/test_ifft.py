from netCDF4 import Dataset
import numpy as np
import modules
import matplotlib.pyplot as plt
from scipy import signal
import numpy.fft as FFT
import scipy.signal
from math import radians, cos, sin, asin, sqrt 
from scipy.fftpack import fft,fftfreq, rfft

#load in data here
#####
#####

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
    
def remove_periodic(bp_periods,bp_mag,ofac):
    
    rm_periods = [1./4,1./3,1./2,1.,365.25/4,365.25/3,365.25/2,365.25]
    rm_n_points = [200,150,100,50,4,3,1,1]
        
    rm_points = []
    
    for i in range(len(rm_periods)):
        rm_i = min(range(len(bp_periods)), key=lambda x: abs(bp_periods[x]-rm_periods[i]))
        n_points = int((rm_n_points[i]*ofac)+np.floor(ofac/2.))
        rm_points_solo = range(rm_i-n_points,(rm_i+n_points)+1)
        if rm_points_solo[0] < 0:
            rm_points_solo = range(0,(rm_i+n_points)+1)
        rm_points = np.append(rm_points,rm_points_solo)
    
    #can either interpolate spectra using points either side of cut, set points to be 0 or delete
    #if going to delete points need to delete equivalent periods too
    
    rm_points = rm_points.astype(int)
    
    #1.
    #bp_mag[rm_points] = np.average((bp_mag[rm_points[0]-1],bp_mag[rm_points[-1]+1]))
    #2.
    bp_mag[rm_points] = 0
    #3.
    #bp_mag = np.delete(bp_mag,rm_points)
    #bp_periods = np.delete(bp_periods,rm_points)
    
    return bp_periods,bp_mag
    
def make_fftarray_from_realimag(fr,fi,data,ofac):
    
    fft_array = [0]*((len(fr)*2)+2) 
    #set first real component to average, imag to 0 
    fft_array[0] = complex((np.average(data)*len(data)*ofac),0)
    
    #Get reverse real and imaginary values
    rev_fr=np.copy(fr[::-1])
    rev_fi=np.copy(fi[::-1])

    f_index = 1

    #Fill half of Fourier Spectrum with positive real and imaginary values
    for i in range(len(fr)):
        fft_array[f_index] = complex(fr[i],fi[i])
        f_index+=1

    fft_array[f_index] = complex(0,0)
    f_index+=1
    
    #Fill Other half of Fourier Spectrum with negative real and imaginary values
    for i in range(len(fr)):
        fft_array[f_index] = complex(rev_fr[i],-rev_fi[i])
        f_index+=1

    #convert to numpy array
    fft_array = np.array(fft_array) 

    return fft_array

def fft_ifft(time,raw_ts,fft_array,ofac,w=False):  
   
    #Take ifft and just take real values
    ifft_ts = np.fft.ifft(fft_array)
    ifft_ts = ifft_ts.astype('float64')

    #cut reconstructed time series if oversampled
    if ofac > 1:
        actual_len = len(time)
        ifft_ts = ifft_ts[:actual_len]
        
        if w == True:
            window = signal.hanning(len(raw_ts))
            ifft_ts = ifft_ts/window
            #cut off 2 percent of ends if windowed.
            len_cut = int((len(ifft_ts)/100.)*2.)
            ifft_ts = ifft_ts[len_cut:-len_cut]
            time = time[len_cut:-len_cut]
            ifft_ts = ifft_ts + np.average(raw_ts)
    
    return time,ifft_ts
    
#take fft
ofac = 1
samp_r = 1./24
fft_periods,fft_mag,fft_ph,fft_fr,fft_fi,fft_array,amp_corr = take_fft(time,data,ofac,samp_r)

#alter fft real/imag components
fft_periods,fft_fr = remove_periodic(fft_periods,fft_fr,ofac)
fft_periods,fft_fi = remove_periodic(fft_periods,fft_fi,ofac)

#remake fft_array from altered fft components
fft_array = make_fftarray_from_realimag(fft_fr,fft_fi,data,ofac)

#take ifft to get time series
ifft_time, ifft_ts = fft_ifft(time,data,fft_array,ofac)

#plot
plt.plot(time,data,color='black')
plt.plot(ifft_time,ifft_ts,color='red')
plt.show()


