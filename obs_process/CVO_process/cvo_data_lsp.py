import numpy as np
import modules
import numpy.fft
import matplotlib.pyplot as plt
import scipy.signal
import redfit
import glob
import csv
from scipy.fftpack import fft,fftfreq, rfft
from scipy import stats

#load in hydrocarbon data

data = np.load('CVO_ALLSPECIES_HOURLY_2006_2014.npy')

time = data[:,2].astype('float64')
ethane = data[:,3].astype('float64')
propane = data[:,4].astype('float64')
isobutane = data[:,5].astype('float64')
nbutane = data[:,6].astype('float64')
acetylene = data[:,7].astype('float64')
isopentane = data[:,8].astype('float64')
npentane = data[:,9].astype('float64')

ethane_orig = np.copy(propane)


gap_test = propane >= 0
cut_time = time[gap_test]
cut_ethane = propane[gap_test]

test = propane < 0
propane[test] = np.nan
not_nan = np.logical_not(np.isnan(propane))
indices = np.arange(len(propane))
ethane = np.interp(indices, indices[not_nan], propane[not_nan])
ethane_ave = np.average(propane)

plt.plot(time,ethane)
plt.show()

#test = propane < 0
#propane[test] = np.nan
#not_nan = np.logical_not(np.isnan(propane))
#indices = np.arange(len(propane))
#propane = np.interp(indices, indices[not_nan], propane[not_nan])
#propane_ave = np.average(propane)

#ethane_time = time[test]
#ethane = ethane[test]
#ethane_ave = np.average(ethane)

#test = propane < 0
#propane_time = time[test]
#propane = propane[test]
#propane_ave = np.average(propane)
#-------------------------------------------

ofac = 1

#ethane_periods,ethane_mag,ethane_ph,ethane_fr,ethane_fi = modules.take_lomb_unwindowed(ethane_time,ethane,ofac)

#closest_ha_index = min(range(len(ethane_periods)), key=lambda i: abs(ethane_periods[i]-182.625))
#closest_annual_index = min(range(len(ethane_periods)), key=lambda i: abs(ethane_periods[i]-365.25))

#rm_indices = [closest_ha_index,closest_annual_index]

#alt_ethane_mag,alt_ethane_fr,alt_ethane_fi,ethane_rm_inds = redfit.sidelobe_n_remove(np.copy(ethane_mag),ethane_fr,ethane_fi,rm_indices,5,ethane_periods)

#plt.loglog(ethane_periods,ethane_mag)
#plt.loglog(ethane_periods,alt_ethane_mag)
#plt.show()

#alt_ethane_mag = np.copy(ethane_mag)
#alt_ethane_fr = np.copy(ethane_fr)
#alt_ethane_fi = np.copy(ethane_fi)

#plt.loglog(ethane_periods,ethane_mag)
#plt.show()


#------------------------------------------------------------------
# #Do IFFT of altered spectra - with significant periods removed and gaps left in real and imag components linearly interpolated.
# #altered spectra provides red noise estimation baseline

# ##use ifft to get time series back from adjusted spectra1
# #complex Fourier spectrum which corresponds to the Lomb-Scargle periodogram:

SAMP_R = 1./24
fft_periods,fft_mag,fft_fr,fft_fi,F2 = modules.take_fft_unwindowed(time,ethane,SAMP_R)

closest_ha_index = min(range(len(fft_periods)), key=lambda i: abs(fft_periods[i]-182.625))
closest_annual_index = min(range(len(fft_periods)), key=lambda i: abs(fft_periods[i]-365.25))
rm_indices = [closest_ha_index,closest_annual_index]

alt_ethane_mag,alt_ethane_fr,alt_ethane_fi,ethane_rm_inds = redfit.sidelobe_n_remove(np.copy(fft_mag),fft_fr,fft_fi,rm_indices,1,fft_periods)


#alt_ethane_mag = np.copy(fft_mag)
#alt_ethane_fr = np.copy(fft_fr)
#alt_ethane_fi = np.copy(fft_fi)

plt.loglog(fft_periods,alt_ethane_mag)
plt.show()


F = [0]*len(F2)

print len(fft_fr)
print len(F)
# 
# # #set first real value to average 
#F[0] = complex(ethane_ave*len(ethane),0)
# 
# # #Get reverse real and imaginary values
#rev_alt_ethane_fr=np.copy(alt_ethane_fr[::-1])
#rev_alt_ethane_fi=np.copy(alt_ethane_fi[::-1])
# 
# 
f_index = 0
# 
# #Fill Fourier Spectrum real and imaginary values
for i in range(len(alt_ethane_fr)):
    F[f_index] = complex(alt_ethane_fr[i],alt_ethane_fi[i])
    f_index+=1
# 
#for i in range(len(alt_ethane_fr)):
#    F[f_index] = complex(rev_alt_ethane_fr[i],-rev_alt_ethane_fi[i])
#    f_index+=1
 
F = np.array(F)    

ifft_ethane = numpy.fft.ifft(F)
ifft_ethane = ifft_ethane.astype('float64')
time = time[gap_test]
ifft_ethane = ifft_ethane[gap_test]

slope, intercept, r_value, p_value, std_err = stats.linregress(time,ifft_ethane)
print intercept,slope
line = slope*time+intercept

plt.plot(cut_time,cut_ethane,linestyle='None',marker='x')
plt.plot(time,ifft_ethane,linestyle='None',marker='x')
plt.plot(time,line)
#plt.plot(time,ifft_ethane)

#all_data = np.vstack((time,ethane_orig,ifft_ethane))
#all_data = np.transpose(all_data)

#speciesheader = ['DAYS_SINCE_2006','ETHANE ORIGINAL(ppbv)','ETHANE DESEASONAL(ppbv)']
#b = open('ethane_deseasonal.txt', 'w')
#a = csv.writer(b)
#a.writerow(speciesheader)
#a.writerows(all_data)
#b.close()


#ifft_ethane_len = (len(ifft_ethane)/ofac) + np.mod(len(ifft_ethane),ofac)

#ifft_time = ethane_time[-ifft_ethane_len:]
#ifft_ethane = ifft_ethane[-len(ethane):]

#ifft_time = ethane_time[-ifft_ethane_len:]
#ifft_ethane = ifft_ethane[-len(ethane):]

#ifft_time = ifft_time-ifft_time[0]

plt.show()

