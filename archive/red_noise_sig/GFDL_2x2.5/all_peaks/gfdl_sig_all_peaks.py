import redfit
import numpy as np
import matplotlib.pyplot as plt
import modules
import scipy.stats
import scipy.signal
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
import numpy.fft
import glob
import sys
import multiprocessing
import operator
import collections
import re

def sort_nicely( l ):
    """ Sort the given list in the way that humans expect."""
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )
    return l

#sort list of GFDL gridbox files into correct order
files = glob.glob('gfdl_sig_all_peaks/grid_files_monthly/*')
#files = ['gfdl_sig_test_2/grid_files_monthly/1817.npy','gfdl_sig_test_2/grid_files_monthly/1819.npy','gfdl_sig_test_2/grid_files_monthly/1860.npy','gfdl_sig_test_2/grid_files_monthly/1890.npy','gfdl_sig_test_2/grid_files_monthly/1811.npy']
#files = ['grid_files/lat-7.0lon181.25.npy','grid_files/lat-7.0lon353.75.npy','grid_files/lat-7.0lon176.25.npy','grid_files/lat-45.0lon188.75.npy']

files = sort_nicely(files)
print files[0:100]

#create significance dictionaries
dict_80_count = {}
dict_85_count = {}
dict_90_count = {}
dict_95_count = {}
dict_99_count = {}
dict_99_9_count = {}

#sorted_files = []
#remove_chars = ['lat','.npy','grid_files/']
#for filename in files:
#    for char in remove_chars:
#        try:
#            stripped = stripped.replace(char,"")
#        except:
#            stripped = filename.replace(char,"")

#    sorted_files.append(stripped)
#    del stripped

#count = 0 
#big_lat_array = np.zeros(len(files))
#big_lon_array = np.zeros(len(files))
#count_array = range(len(files))

#for filename in sorted_files:
#    latlon_array = filename.split('lon')
#    big_lat_array[count] = float(latlon_array[0])
#    big_lon_array[count] = float(latlon_array[1])
#    count+=1
    
#big_latlon_array = np.vstack((big_lat_array,big_lon_array,count_array))
#names = ['lats','lons','inds']
#big_latlon_array = np.core.records.fromarrays(big_latlon_array, names='lats, lons, inds', formats = 'f8, f8, i8')
#big_latlon_array = np.sort(big_latlon_array, order=['lats', 'lons']) 
#sorted_indices = big_latlon_array['inds']
#files = np.array(files)
#files = files[sorted_indices]

print len(files)
grid_count = 0
def read(files,j): 
    print j
    #load in gridbox
    y = np.load(files[j])
    #make ppbv
    #x = x*1e9

    #daily average data
    #daily_x = np.zeros(len(x)/24.) 
    #for i in range(len(daily_x)):
    #    start = int(24*(i))
    #    end = int(24*(i+1))
    #    daily_x[i] = np.average(x[start:end])

    #monthly average data
    #val =365./12

    #season_x = np.zeros(len(daily_x)/(val)) 
    #for i in range(len(season_x)):
    #    start = int(val*(i))
    #    end = int(val*(i+1))
    #    season_x[i] = np.average(daily_x[start:end])

    #y = season_x
    #np.save('gfdl_sig_test/grid_files_monthly/%i'%(j),y)
    #x = np.arange(len(season_x))
    #x = x*val #t in days
    
    #try:
    #    all_y = np.vstack((all_y,y))
    #except:
    #    all_y = y
    return y

if __name__ == '__main__':
    pool = multiprocessing.Pool(processes=32)
    results = [pool.apply_async(read, (files,i)) for i in range(len(files))]
    results = [r.get() for r in results]
    all_y = np.array(results)

print all_y.shape
val =365./12        
x = np.arange(len(all_y[0,:]))
x = x*val #t in days
def main_arg(x,y,grid_count):  
    #set key parameters
    #mctest = [True/False] Calc. MC-based false-alarm levels (Default False)
    #nsim = Number of simulations
    #n50 = number of WOSA segments (50 % overlap)

    mctest = True
    nsim = 50
    n50 = 3
    ofac = 4

    print x
    print y
    #synthetic data
    #pi2 = np.pi*2
    #x = np.arange(0,1001,1)
    #y = 50 + (10*(np.sin((pi2*x/10)-(0))))

    var_mean = np.mean(y)
    #subtract mean from data
    y = y - var_mean

    #average dt of entire time series
    diffs = [x[i+1]-x[i] for i in range(len(x)-1)]  
    avgdt = np.average(diffs)

    #take Lomb-Scargle
    var_mean = np.mean(y)
    #periods,freqs,mag,full_n,fr,fi,ave_seg_len = redfit.welch_lomb(x,y,ofac,n50,avgdt,detrend=False)
    periods,mag,ph,sig,i_freq,fr,fi = modules.take_lomb(x,y,ofac,avgdt)
    freqs = 1./periods

    #find indices of strong periodicities					    
    closest_annual = min(range(len(periods)), key=lambda i: abs(periods[i]-365))
    closest_hannual = min(range(len(periods)), key=lambda i: abs(periods[i]-365./2))
    closest_third = min(range(len(periods)), key=lambda i: abs(periods[i]-365./3))
    closest_quarter = min(range(len(periods)), key=lambda i: abs(periods[i]-365./4))
    closest_daily = min(range(len(periods)), key=lambda i: abs(periods[i]-1))
	
    #remove significant periods

    annual_test = (periods >= 350) & (periods <= 380)
    hannual_test = (periods >= 180) & (periods <= 185)
    third_test = (periods >= 120) & (periods <= 123)
    quarter_test = (periods >= 90) & (periods <= 92)

    corr_mag = np.copy(mag)
    corr_fr = np.copy(fr)
    corr_fi = np.copy(fi)

    set_annual_i = (np.argmax(annual_test==True))-1
    set_hannual_i = (np.argmax(hannual_test==True))-1
    set_third_i = (np.argmax(third_test==True))-1
    set_quarter_i = (np.argmax(quarter_test==True))-1

    set_annual_i = np.arange(set_annual_i-len(corr_mag[annual_test]),set_annual_i)
    set_hannual_i = np.arange(set_hannual_i-len(corr_mag[hannual_test]),set_hannual_i)
    set_third_i = np.arange(set_third_i-len(corr_mag[third_test]),set_third_i)
    set_quarter_i = np.arange(set_quarter_i-len(corr_mag[quarter_test]),set_quarter_i)

    annual_mag_set_val = corr_mag[set_annual_i]
    hannual_mag_set_val = corr_mag[set_hannual_i]
    third_mag_set_val = corr_mag[set_third_i]
    quarter_mag_set_val = corr_mag[set_quarter_i]

    annual_fr= corr_fr[set_annual_i]
    annual_fi= corr_fi[set_annual_i]
    hannual_fr= corr_fr[set_hannual_i]
    hannual_fi= corr_fi[set_hannual_i]
    third_fr= corr_fr[set_third_i]
    third_fi= corr_fi[set_third_i]
    quarter_fr= corr_fr[set_quarter_i]
    quarter_fi= corr_fi[set_quarter_i]

    corr_mag[annual_test] = annual_mag_set_val
    corr_fr[annual_test] = annual_fr
    corr_fi[annual_test] = annual_fi

    corr_mag[hannual_test] = hannual_mag_set_val
    corr_fr[hannual_test] = hannual_fr
    corr_fi[hannual_test] = hannual_fr

    corr_mag[third_test] = third_mag_set_val
    corr_fr[third_test] = third_fr
    corr_fi[third_test] = third_fi

    corr_mag[quarter_test] = quarter_mag_set_val
    corr_fr[quarter_test] = quarter_fr
    corr_fi[quarter_test] = quarter_fi

#------------------------------------------------------
    #Section uses ifft to generate time series of detrended, key period treated data
    # make complex Fourier spectrum which corresponds to the Lomb-Scargle periodogram:
    F = [0]*((len(corr_fr)*2)+1)
    #set first real value to average 
    F[0] = complex(var_mean,0)

    #Get reverse real and imaginary values
    rev_fr=corr_fr[::-1] 
    rev_fi=corr_fi[::-1]

    #Fill Fourier Spectrum real and imaginary values
    for i in range(len(corr_fr)):
        F[i+1] = complex(corr_fr[i],corr_fi[i])
    for i in range(len(corr_fr),len(corr_fr)*2):
        F[i+1] = complex(rev_fr[i-len(corr_fr)],-rev_fi[i-len(corr_fr)])
    F = np.array(F)

    #Take ifft and just take real values
    output_ts = numpy.fft.ifft(F)
    output_ts = output_ts.astype('float64')

    time = []

    #generate time points
    #	for i in range(len(F)):
    #    time= np.append(time,x[0] +(i/(freqs[1]-freqs[0])/len(F)))

    time = x

    time = time[-len(y):]
    output_ts = output_ts[-len(y):]

    #gap = int((len(output_ts) - ave_seg_len)/2)
    #output_ts = output_ts[
    output_ts = output_ts/np.hanning(len(output_ts))
    time = time[350:-350]
    output_ts = output_ts[350:-350]
    #time = time[gap:-gap]
    #output_ts = output_ts[gap:-gap]
    output_ts = output_ts+var_mean

    #save out modified time series
    np.save('gfdl_sig_all_peaks/modified_GFDL_ts/%i'%(grid_count),output_ts)
    
    #------------------------------------------------------
    #Section uses ifft to generate time series of detrended, key period treated data
    # make complex Fourier spectrum which corresponds to the Lomb-Scargle periodogram:
    F = [0]*((len(fr)*2)+1)
    #set first real value to average 
    F[0] = complex(var_mean,0)

    #Get reverse real and imaginary values
    rev_fr=fr[::-1] 
    rev_fi=fi[::-1]

    #Fill Fourier Spectrum real and imaginary values
    for i in range(len(fr)):
        F[i+1] = complex(fr[i],fi[i])
    for i in range(len(fr),len(fr)*2):
        F[i+1] = complex(rev_fr[i-len(fr)],-rev_fi[i-len(fr)])
    F = np.array(F)

    #Take ifft and just take real values
    other_output_ts = numpy.fft.ifft(F)
    other_output_ts = other_output_ts.astype('float64')

    other_time = []

    #generate time points
    #	for i in range(len(F)):
    #    time= np.append(time,x[0] +(i/(freqs[1]-freqs[0])/len(F)))

    other_time = x

    other_time = other_time[-len(y):]
    other_output_ts = other_output_ts[-len(y):]

    #gap = int((len(output_ts) - ave_seg_len)/2)
    #output_ts = output_ts[
    other_output_ts = other_output_ts/np.hanning(len(other_output_ts))
    other_time = other_time[350:-350]
    other_output_ts = other_output_ts[350:-350]
    #time = time[gap:-gap]
    #output_ts = output_ts[gap:-gap]
    other_output_ts = other_output_ts+var_mean
    model_periods,model_freq,model_mag,full_n,fr,fi,ave_seg_len = redfit.welch_lomb(other_time,other_output_ts,ofac,n50,avgdt,detrend=True)

#-------------------------------------------------------
    alt_periods,alt_mag,red_periods,red_mag,gredth,fac80,fac85,fac90,fac95,fac99,fac99_9,faccrit,tau,corr = redfit.redcore(nsim,n50,mctest,time,output_ts)
    
    #correct unmodified spectrum
    model_mag = model_mag/corr
    if tau == 'invalid':
        print 'Tau Invalid.'
        np.save('gfdl_sig_all_peaks/invalid_gridboxes/%i'%(grid_count),output_ts)
        mag_80,mag_85,mag_90,mag_95,mag_99,mag_99_9,periods_80,periods_85,periods_90,periods_95,periods_99,periods_99_9 = [],[],[],[],[],[],[],[],[],[],[],[]
        return mag_80,mag_85,mag_90,mag_95,mag_99,mag_99_9,periods_80,periods_85,periods_90,periods_95,periods_99,periods_99_9
    
    #save out output spectrums
    np.save('gfdl_sig_all_peaks/modified_GFDL_spectra/%i'%(grid_count),model_mag)
    np.save('gfdl_sig_all_peaks/periods',model_periods)

    #fig=plt.figure(figsize=(16,8))
    #fig.patch.set_facecolor('white')
    #ax = fig.add_subplot(1,1,1)

    #def form2(x, pos):
    #	""" This function returns a string with 3 decimal places, given the input x"""
    #	return '%.2f' % x	

    #def form5(x, pos):
    #	""" This function returns a string with 3 decimal places, given the input x"""
    #	return '%.5f' % x

    #ax.grid(True)
    #plt.loglog(model_periods,model_mag,color='black', label = 'Obs. Spectrum')
    #plt.loglog(red_periods,red_mag,color='red', label = 'Red Noise')
    #plt.loglog(red_periods,gredth,color='red', label = 'Theoretical Red Noise Spectrum')
    #plt.loglog(red_periods,gredth*fac80,color='orange', label = '80% Significance')
    #plt.loglog(red_periods,gredth*fac85,color='blue',label = '85% Significance')
    #plt.loglog(red_periods,gredth*fac90,color='purple', label = '90% Significance')
    #plt.loglog(red_periods,gredth*fac95,color='orange', label = '95% Significance')
    #plt.loglog(red_periods,gredth*fac99,color='blue',label = '99% Significance')
    #plt.loglog(red_periods,gredth*fac99_9,color='purple', label = '99.9% Significance')
    #leg=ax.legend(loc=0, prop={'size':15})
    #leg.get_frame().set_alpha(0.4)

    #xformatter = FuncFormatter(form2)
    #yformatter = FuncFormatter(form2)
    #ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
    #plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
    #plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))

    #plt.show()

#------------------------
    #Determine if spectal values are significant at various significance levels

    test_80 = model_mag >= gredth*fac80
    test_85 = model_mag >= gredth*fac85
    test_90 = model_mag >= gredth*fac90
    test_95 = model_mag >= gredth*fac95
    test_99 = model_mag >= gredth*fac99
    test_99_9 = model_mag >= gredth*fac99_9

    periods_80 = model_periods[test_80]
    periods_85 = model_periods[test_85]
    periods_90 = model_periods[test_90]
    periods_95 = model_periods[test_95]
    periods_99 = model_periods[test_99]
    periods_99_9 = model_periods[test_99_9]

    mag_80 = model_mag[test_80]
    mag_85 = model_mag[test_85]
    mag_90 = model_mag[test_90]
    mag_95 = model_mag[test_95]
    mag_99 = model_mag[test_99]
    mag_99_9 = model_mag[test_99_9]
        
    grid_count+=1
    
    return mag_80,mag_85,mag_90,mag_95,mag_99,mag_99_9,periods_80,periods_85,periods_90,periods_95,periods_99,periods_99_9

all_sig_mag_80 = [0]*len(files)
all_sig_mag_85 = [0]*len(files)
all_sig_mag_90 = [0]*len(files)
all_sig_mag_95 = [0]*len(files)
all_sig_mag_99 = [0]*len(files)
all_sig_mag_99_9 = [0]*len(files)
all_sig_period_80 = [0]*len(files)
all_sig_period_85 = [0]*len(files)
all_sig_period_90 = [0]*len(files)
all_sig_period_95 = [0]*len(files)
all_sig_period_99 = [0]*len(files)
all_sig_period_99_9 = [0]*len(files)

chunk = len(files)/5
for j in range(5):
    if __name__ == '__main__':
        chunk_start = chunk*j
        chunk_end = chunk*(j+1)
        pool = multiprocessing.Pool(processes=32)
        results = [pool.apply_async(main_arg, (x,all_y[i],i)) for i in range(chunk_start,chunk_end)]
        results = [r.get() for r in results]
        results = np.array(results)
    
        sig_mag_80 = results[:,0]
        sig_mag_85 = results[:,1]
        sig_mag_90 = results[:,2]
        sig_mag_95 = results[:,3]
        sig_mag_99 = results[:,4]
        sig_mag_99_9 = results[:,5]
        sig_period_80 = results[:,6]
        sig_period_85 = results[:,7]
        sig_period_90 = results[:,8]
        sig_period_95 = results[:,9]
        sig_period_99 = results[:,10]
        sig_period_99_9 = results[:,11]
        
        all_sig_mag_80[chunk_start:chunk_end] = sig_mag_80
        all_sig_mag_85[chunk_start:chunk_end] = sig_mag_85
        all_sig_mag_90[chunk_start:chunk_end] = sig_mag_90
        all_sig_mag_95[chunk_start:chunk_end] = sig_mag_95
        all_sig_mag_99[chunk_start:chunk_end] = sig_mag_99
        all_sig_mag_99_9[chunk_start:chunk_end] = sig_mag_99_9
        
        all_sig_period_80[chunk_start:chunk_end] = sig_period_80
        all_sig_period_85[chunk_start:chunk_end] = sig_period_85
        all_sig_period_90[chunk_start:chunk_end] = sig_period_90
        all_sig_period_95[chunk_start:chunk_end] = sig_period_95
        all_sig_period_99[chunk_start:chunk_end] = sig_period_99
        all_sig_period_99_9[chunk_start:chunk_end] = sig_period_99_9

flat_sig_period_80 = []
flat_sig_period_85 = []
flat_sig_period_90 = []
flat_sig_period_95 = []
flat_sig_period_99 = []
flat_sig_period_99_9 = []
        
for i in all_sig_period_80:
    flat_sig_period_80 = np.append(flat_sig_period_80,i)
for i in all_sig_period_85:
    flat_sig_period_85 = np.append(flat_sig_period_85,i)
for i in all_sig_period_90:
    flat_sig_period_90 = np.append(flat_sig_period_90,i)
for i in all_sig_period_95:
    flat_sig_period_95 = np.append(flat_sig_period_95,i)
for i in all_sig_period_99:
    flat_sig_period_99 = np.append(flat_sig_period_99,i)
for i in all_sig_period_99_9:
    flat_sig_period_99_9 = np.append(flat_sig_period_99_9,i)
    
count_period_80 =  collections.Counter(flat_sig_period_80)    
count_period_85 =  collections.Counter(flat_sig_period_85)
count_period_90 =  collections.Counter(flat_sig_period_90)
count_period_95 =  collections.Counter(flat_sig_period_95)
count_period_99 =  collections.Counter(flat_sig_period_99)
count_period_99_9 =  collections.Counter(flat_sig_period_99_9)

np.save('gfdl_sig_all_peaks/significant_mags/sig_mag_80',all_sig_mag_80)
np.save('gfdl_sig_all_peaks/significant_mags/sig_mag_85',all_sig_mag_85)
np.save('gfdl_sig_all_peaks/significant_mags/sig_mag_90',all_sig_mag_90)
np.save('gfdl_sig_all_peaks/significant_mags/sig_mag_95',all_sig_mag_95)
np.save('gfdl_sig_all_peaks/significant_mags/sig_mag_99',all_sig_mag_99)
np.save('gfdl_sig_all_peaks/significant_mags/sig_mag_99_9',all_sig_mag_99_9)
    
np.save('gfdl_sig_all_peaks/significant_periods/sig_period_80',all_sig_period_80)
np.save('gfdl_sig_all_peaks/significant_periods/sig_period_85',all_sig_period_85)
np.save('gfdl_sig_all_peaks/significant_periods/sig_period_90',all_sig_period_90)
np.save('gfdl_sig_all_peaks/significant_periods/sig_period_95',all_sig_period_95)
np.save('gfdl_sig_all_peaks/significant_periods/sig_period_99',all_sig_period_99)
np.save('gfdl_sig_all_peaks/significant_periods/sig_period_99_9',all_sig_period_99_9)
    
np.save('gfdl_sig_all_peaks/counts/period_80',count_period_80)
np.save('gfdl_sig_all_peaks/counts/period_85',count_period_85)
np.save('gfdl_sig_all_peaks/counts/period_90',count_period_90)
np.save('gfdl_sig_all_peaks/counts/period_95',count_period_95)
np.save('gfdl_sig_all_peaks/counts/period_99',count_period_99)
np.save('gfdl_sig_all_peaks/counts/period_99_9',count_period_99_9)
 
