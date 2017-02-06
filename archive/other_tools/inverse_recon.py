import matplotlib.pyplot as plt 
from scipy.fftpack import fftfreq,fft, ifft
import numpy as np
from cmath import *
import operator
import lomb_phase
import lomb
import modules

#Time array in hourly intervals

hourly_step = 1./24

a = np.arange(0,365.25*4,hourly_step)


pi2 = np.pi*2

daily_amplitude = 2
half_annual_amplitude = 5
annual_amplitude = 10

cos_waveform1 = daily_amplitude*(np.cos((pi2*a/1)-0))
cos_waveform2 = half_annual_amplitude*(np.cos((pi2*a/(365.25/2))-0))
cos_waveform3 = annual_amplitude*(np.sin((pi2*a/365.25)+(0)))

waveform = 50 + (cos_waveform1+cos_waveform2+cos_waveform3)
#waveform = 50 + (cos_waveform3)

#do FFT

fft_array = fft(waveform)
#fft_array = fft_array/(1./hourly_step)
fft_real = fft_array.real
fft_imag = fft_array.imag
#fft_cos_amp = np.array(fft_cos_amp)
#fft_sin_amp = np.array(fft_sin_amp)

fft_mag = np.abs(fft_array)
fft_phase = [phase(i) for i in fft_array]
fft_phase = np.array(fft_phase)
fft_freqs = fftfreq(len(a),hourly_step) 

#valid = fft_freqs > 0
#fft_mag, fft_phase, fft_freqs = fft_mag[valid], fft_phase[valid], fft_freqs[valid] 

#fft_real = fft_imag/(len(a))
#fft_imag = fft_imag/(len(a))
#fft_mag = fft_mag/(len(a)/2)
#fft_periods = 1./fft_freqs

output = ifft( fft( waveform ) )

#Process fft array

pos_daily_index = np.argmin(np.abs(fft_freqs - 1))
pos_half_annual_index = np.argmin(np.abs(fft_freqs - 1./(365.25/2)))
pos_annual_index = np.argmin(np.abs(fft_freqs - (1/365.25)))

neg_daily_index = np.argmin(np.abs(fft_freqs - (-1)))
neg_half_annual_index = np.argmin(np.abs(fft_freqs - 1./(-365.25/2)))
neg_annual_index = np.argmin(np.abs(fft_freqs - (1/-365.25)))


real_indices = [0,pos_daily_index,neg_daily_index,pos_half_annual_index,neg_half_annual_index,pos_annual_index,neg_annual_index]
imag_indices=[pos_daily_index,neg_daily_index,pos_half_annual_index,neg_half_annual_index,pos_annual_index,neg_annual_index]


mask = np.ones(fft_real.shape,dtype=bool) #np.ones_like(a,dtype=bool)
mask[real_indices] = False
fft_real[mask] = 0

mask = np.ones(fft_imag.shape,dtype=bool) #np.ones_like(a,dtype=bool)
mask[imag_indices] = False
fft_imag[mask] = 0


edited_fft_array = np.vectorize(complex)(fft_real,fft_imag)


inverse_fft =  ifft(edited_fft_array) 


#plot up

fig = plt.figure(figsize=(14, 10 ))
fig.patch.set_facecolor('white')
ax1 = plt.subplot(211)
#plt.plot(a,cos_waveform1)
#plt.plot(a,cos_waveform2)
#plt.plot(a,cos_waveform3)
#ax1.plot(a,waveform)
#ax1.semilogx(fft_periods, fft_mag)

ax2 = plt.subplot(212)
#ax2.semilogx(fft_periods, fft_mag)
#ax2.set_xlim([0.9,5])
#ax2.set_ylim([0, 10])
#plt.grid()
#plt.xlabel('Period (Days)', fontsize = 16)
#plt.ylabel('Magnitude (ppbV)', fontsize = 16)
#plt.title('Lomb-Scargle derived FFT', fontsize = 16)
#plt.plot(a, waveform)
#plt.plot(data_points,data)
#ax2.plot(fft_periods,fft_mag)
ax2.plot(a, inverse_fft )
#plt.ylim(40,60)
#plt.plot(combined_time, waveform)
#plt.semilogx(fft_periods, fft_sin_amp)
#plt.semilogx(full_freqs, full_real)
#plt.loglog(fft_periods,fft_mag)
#plt.semilogx(lomb_periods,lomb_mag)
#plt.ylim(0,15)

#plt.axvline(1, color = 'red', linestyle = '--')
#plt.axvline(365.25/2, color = 'red', linestyle = '--')
#plt.axvline(365.25, color = 'red', linestyle = '--')

#plt.text(0.6,14,"Daily Period")
#plt.text(50,13,"Half Annual Period")
#plt.text(230,14,"Annual Period")

plt.show()
