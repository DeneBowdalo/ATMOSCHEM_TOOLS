import matplotlib.pyplot as plt
from netCDF4 import Dataset
import numpy as np
import matplotlib
import modules

kind = raw_input('\namp or ph?\n')
period = 's'

#raw_obs_data = Dataset('/work/home/db876/observations/ozonesonde/process/O3_RADIOSONDES_2005_2010.nc')
obs_data = Dataset('obs_sig_periods.nc')
model_data = Dataset('../GEOSCHEM_VERT_v902_2x2.5_GEOS5_D_*/model_sig_periods.nc')

lats = obs_data.variables['latitude'][:]
lons = obs_data.variables['longitude'][:]
countries = obs_data.variables['countries'][:]
refs = obs_data.variables['refs'][:]
pressures = obs_data.variables['pressure_centres'][:]

obs_seasonal_waveforms  = obs_data.variables['seasonal_waveform'][:]
obs_seasonal_amplitudes  = obs_data.variables['seasonal_amplitude'][:]
obs_seasonal_phases  = obs_data.variables['seasonal_max_phase'][:]
model_seasonal_waveforms = model_data.variables['seasonal_waveform'][:]
model_seasonal_amplitudes  = model_data.variables['seasonal_amplitude'][:]
model_seasonal_phases  = model_data.variables['seasonal_max_phase'][:] 

test = obs_seasonal_waveforms < 0
obs_seasonal_waveforms[test] = np.NaN

test = obs_seasonal_amplitudes < 0
obs_seasonal_amplitudes[test] = np.NaN

test = model_seasonal_waveforms < 0
model_seasonal_waveforms[test] = np.NaN

area_boundaries,area_tags,area_labels = modules.area_dicts()
tags = modules.get_tags(refs)

areas = ['ANT','S_O','AF','N_O','S_EU','C_EU','NE_NA','ARC']

obsn_linspace = np.linspace(0, 1, len(areas))

fig, axes = plt.subplots(nrows=2, ncols=4,figsize=(19,13))
fig.patch.set_facecolor('white')
count = 0
for ax in axes.flat:
    try:
        area = areas[count]
    except:
        ax.axis('off')
        continue

    area = areas[count]

    area_grid = area_boundaries[area]
    area_tag = area_tags[area]
    area_label = area_labels[area]
    cut_test = modules.area_cut(area,lats,lons,tags,area_grid,area_tag)

    cut_refs = refs[cut_test]
    
    y_array = np.average(pressures[cut_test,:],axis=0)

    print obs_seasonal_waveforms.shape
    obs_seasonal_waveform_cut = obs_seasonal_waveforms[cut_test,:,:]
    obs_seasonal_amp_cut = obs_seasonal_amplitudes[cut_test,:]
    obs_seasonal_ph_cut = obs_seasonal_phases[cut_test,:]
    model_seasonal_waveform_cut = model_seasonal_waveforms[cut_test,:,:]
    model_seasonal_amp_cut = model_seasonal_amplitudes[cut_test,:]
    model_seasonal_ph_cut = model_seasonal_phases[cut_test,:]
    print obs_seasonal_waveform_cut.shape

    print len(cut_refs)
    for i in range(len(cut_refs)):
    #obs_seasonal_waveform = np.nanmean(obs_seasonal_waveforms[cut_test,:,:],axis=0)
    #model_seasonal_waveform = np.nanmean(model_seasonal_waveforms[cut_test,:,:],axis=0)
        obs_seasonal_waveform = obs_seasonal_waveform_cut[i,:,:]
        model_seasonal_waveform = model_seasonal_waveform_cut[i,:,:]

        obs_seasonal_amp = obs_seasonal_amp_cut[i,:]
        model_seasonal_amp = model_seasonal_amp_cut[i,:]
    
        obs_seasonal_ph = obs_seasonal_ph_cut[i,:]
        model_seasonal_ph = model_seasonal_ph_cut[i,:]

        #get new average amplitude and max ph
        #obs_seasonal_amp = (np.max(obs_seasonal_waveform,axis=0)-np.min(obs_seasonal_waveform,axis=0))/2.
        #print obs_seasonal_amp.shape
        #obs_max_ph = []
        #test = np.isnan(obs_seasonal_waveform)
        #for i in range(len(y_array)):
        #    test = np.isnan(obs_seasonal_waveform[:,i])
        #    if False in test:
        #        obs_max_ph.append(np.argmax(obs_seasonal_waveform[:,i]))
        #    else:
        #        obs_max_ph.append(np.nan)

        #model_seasonal_amp = (np.max(model_seasonal_waveform,axis=0)-np.min(model_seasonal_waveform,axis=0))/2.
        #model_max_ph = []                                                                                                                                                          
        #test = np.isnan(model_seasonal_waveform)
        #for i in range(len(y_array)):
       #     test = np.isnan(model_seasonal_waveform[:,i])
       #     if False in test:
       #         model_max_ph.append(np.argmax(model_seasonal_waveform[:,i]))
       #     else:
       #         model_max_ph.append(np.nan)


        if kind == 'amp':
            print obs_seasonal_amp.shape
            print y_array.shape
            pl = ax.scatter(obs_seasonal_amp,y_array,c='black')
            pl2 = ax.scatter(model_seasonal_amp,y_array,c='red')
            ax.set_xlim(0,30)                                                                                                                                                      
            ax.set_ylim(1000,400)
            #ax.set_xlabel('Seasonal Amplitude',fontsize=20)

        elif kind == 'ph':
            pl = ax.scatter(obs_seasonal_ph,y_array,c='black')
            pl2 = ax.scatter(model_seasonal_ph,y_array,c='red')
            ax.set_xlim(0,12)
            ax.set_ylim(1000,400)
            #ax.set_xlabel('Seasonal Phase',fontsize=20)

        print '\n'

    ax.text(0.81, 0.95, area_label, ha='center', va='center',transform=ax.transAxes,fontsize=15)
    count+=1

    

#ax.set_ylabel('Altitude (m)',fontsize=20)    
#cb = fig.colorbar(pl)
#cb.set_label('Ozone Average (ppb)',fontsize=20)

    #plt.title('%s, %s, %s Latitude, %s Longitude'%(ref,country,lat,lon),fontsize = 20)
    
    #plt.savefig('plots/%s_%s/%s.png'%(descriptor,kind,ref))

    #plt.close()
plt.show()
