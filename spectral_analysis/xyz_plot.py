#%matplotlib notebook 
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from ipywidgets import Dropdown
from IPython.display import display
from ipywidgets import *
import modules

#set up writing to ipython notebook
#output_notebook()

#-----------------------------------------------------
#get species from current directory
present_dir = os.getcwd()
obs_fname,species,start_year,end_year,vres,timeres,run_type = modules.get_obs_info(present_dir)
model_fname,species,start_year,end_year,vres,timeres = modules.get_model_info(present_dir)
    
#read in obs_ts_data
obs_refs,obs_raw_time,obs_ref_time,obs_datetime_time,obs_std_var,obs_lats,obs_lons,obs_alt,obs_process_groups,obs_raw_class,obs_anthrome_class,obs_gap_inds = modules.read_obs_all(obs_fname,species,start_year,end_year)
model_raw_time,model_ref_time,model_datetime_time,model_std_var,lat_e,lon_e,lat_c,lon_c,grid_size,gridbox_count = modules.read_model_all(model_fname,species,start_year,end_year)

#get obs lat_lon grid central points
obs_lats_centre, obs_lons_centre, model_indices = modules.obs_model_gridboxes(lat_e,lon_e,obs_lats,obs_lons) 

#get observational location tags 
tags = modules.get_tags(obs_refs)
loc_colors = {'ANT':'magenta','ARC':'purple','AF':'orange','EA':'blue','MA':'blue','EU':'red','OC':'black','O':'yellow','NA':'green','SA':'pink'}
tag_colors = [loc_colors[tag] for tag in tags]

#--------------------------------------------------------
#load in periodic lsp data
obs_periodic_fname = '../obs_%s_%s/LSP_stats.nc'%(vres,timeres)
model_periodic_fname = 'LSP_stats.nc'

(obs_diurnal_mag,obs_diurnal_ph,obs_seasonal_mag,obs_seasonal_ph,obs_mean,obs_p1,obs_p5,obs_p25,obs_p50,obs_p75,obs_p95,obs_p99,obs_diurnal_ave_waveform,obs_diurnal_ave_waveform_extended,obs_diurnal_season_waveform_extended,obs_seasonal_waveform,obs_seasonal_waveform_extended,obs_full_ave_waveform,obs_full_season_waveform,obs_pc_var_daily,obs_pc_var_seasonal,obs_pc_var_full,
obs_diurnal_mag_spring,obs_diurnal_ph_spring,obs_mean_spring,obs_p1_spring,obs_p5_spring,obs_p25_spring,obs_p50_spring,obs_p75_spring,obs_p95_spring,obs_p99_spring,obs_diurnal_waveform_spring,
obs_diurnal_mag_summer,obs_diurnal_ph_summer,obs_mean_summer,obs_p1_summer,obs_p5_summer,obs_p25_summer,obs_p50_summer,obs_p75_summer,obs_p95_summer,obs_p99_summer,obs_diurnal_waveform_summer,
obs_diurnal_mag_autumn,obs_diurnal_ph_autumn,obs_mean_autumn,obs_p1_autumn,obs_p5_autumn,obs_p25_autumn,obs_p50_autumn,obs_p75_autumn,obs_p95_autumn,obs_p99_autumn,obs_diurnal_waveform_autumn,
obs_diurnal_mag_winter,obs_diurnal_ph_winter,obs_mean_winter,obs_p1_winter,obs_p5_winter,obs_p25_winter,obs_p50_winter,obs_p75_winter,obs_p95_winter,obs_p99_winter,obs_diurnal_waveform_winter,
obs_seasonal_mag_day,obs_seasonal_ph_day,obs_mean_day,obs_p1_day,obs_p5_day,obs_p25_day,obs_p50_day,obs_p75_day,obs_p95_day,obs_p99_day,obs_seasonal_waveform_day,
obs_seasonal_mag_night,obs_seasonal_ph_night,obs_mean_night,obs_p1_night,obs_p5_night,obs_p25_night,obs_p50_night,obs_p75_night,obs_p95_night,obs_p99_night,obs_seasonal_waveform_night,
obs_daily_h3_mag,obs_daily_h2_mag,obs_daily_h1_mag,obs_daily_mag,obs_annual_h3_mag,obs_annual_h2_mag,obs_annual_h1_mag,obs_annual_mag,
model_diurnal_mag,model_diurnal_ph,model_seasonal_mag,model_seasonal_ph,model_mean,model_p1,model_p5,model_p25,model_p50,model_p75,model_p95,model_p99,model_diurnal_ave_waveform,model_diurnal_ave_waveform_extended,model_diurnal_season_waveform_extended,model_seasonal_waveform,model_seasonal_waveform_extended,model_full_ave_waveform,model_full_season_waveform,model_pc_var_daily,model_pc_var_seasonal,model_pc_var_full,
model_diurnal_mag_spring,model_diurnal_ph_spring,model_mean_spring,model_p1_spring,model_p5_spring,model_p25_spring,model_p50_spring,model_p75_spring,model_p95_spring,model_p99_spring,model_diurnal_waveform_spring,
model_diurnal_mag_summer,model_diurnal_ph_summer,model_mean_summer,model_p1_summer,model_p5_summer,model_p25_summer,model_p50_summer,model_p75_summer,model_p95_summer,model_p99_summer,model_diurnal_waveform_summer, 
model_diurnal_mag_autumn,model_diurnal_ph_autumn,model_mean_autumn,model_p1_autumn,model_p5_autumn,model_p25_autumn,model_p50_autumn,model_p75_autumn,model_p95_autumn,model_p99_autumn,model_diurnal_waveform_autumn, 
model_diurnal_mag_winter,model_diurnal_ph_winter,model_mean_winter,model_p1_winter,model_p5_winter,model_p25_winter,model_p50_winter,model_p75_winter,model_p95_winter,model_p99_winter,model_diurnal_waveform_winter,
model_seasonal_mag_day,model_seasonal_ph_day,model_mean_day,model_p1_day,model_p5_day,model_p25_day,model_p50_day,model_p75_day,model_p95_day,model_p99_day,model_seasonal_waveform_day,
model_seasonal_mag_night,model_seasonal_ph_night,model_mean_night,model_p1_night,model_p5_night,model_p25_night,model_p50_night,model_p75_night,model_p95_night,model_p99_night,model_seasonal_waveform_night,
model_daily_h3_mag,model_daily_h2_mag,model_daily_h1_mag,model_daily_mag,model_annual_h3_mag,model_annual_h2_mag,model_annual_h1_mag,model_annual_mag) = modules.get_periodic_specific(obs_periodic_fname,model_periodic_fname)

#group data

#xy_params
xy_params = [obs_diurnal_mag,obs_diurnal_ph,obs_seasonal_mag,obs_seasonal_ph,obs_mean,obs_p1,obs_p5,obs_p25,obs_p50,obs_p75,obs_p95,obs_p99,obs_pc_var_daily,obs_pc_var_seasonal,obs_pc_var_full,
obs_diurnal_mag_spring,obs_diurnal_ph_spring,obs_mean_spring,obs_p1_spring,obs_p5_spring,obs_p25_spring,obs_p50_spring,obs_p75_spring,obs_p95_spring,obs_p99_spring,
obs_diurnal_mag_summer,obs_diurnal_ph_summer,obs_mean_summer,obs_p1_summer,obs_p5_summer,obs_p25_summer,obs_p50_summer,obs_p75_summer,obs_p95_summer,obs_p99_summer,
obs_diurnal_mag_autumn,obs_diurnal_ph_autumn,obs_mean_autumn,obs_p1_autumn,obs_p5_autumn,obs_p25_autumn,obs_p50_autumn,obs_p75_autumn,obs_p95_autumn,obs_p99_autumn,
obs_diurnal_mag_winter,obs_diurnal_ph_winter,obs_mean_winter,obs_p1_winter,obs_p5_winter,obs_p25_winter,obs_p50_winter,obs_p75_winter,obs_p95_winter,obs_p99_winter,
obs_seasonal_mag_day,obs_seasonal_ph_day,obs_mean_day,obs_p1_day,obs_p5_day,obs_p25_day,obs_p50_day,obs_p75_day,obs_p95_day,obs_p99_day,
obs_seasonal_mag_night,obs_seasonal_ph_night,obs_mean_night,obs_p1_night,obs_p5_night,obs_p25_night,obs_p50_night,obs_p75_night,obs_p95_night,obs_p99_night,
obs_daily_h3_mag,obs_daily_h2_mag,obs_daily_h1_mag,obs_daily_mag,obs_annual_h3_mag,obs_annual_h2_mag,obs_annual_h1_mag,obs_annual_mag,
model_diurnal_mag,model_diurnal_ph,model_seasonal_mag,model_seasonal_ph,model_mean,model_p1,model_p5,model_p25,model_p50,model_p75,model_p95,model_p99,model_pc_var_daily,model_pc_var_seasonal,model_pc_var_full,
model_diurnal_mag_spring,model_diurnal_ph_spring,model_mean_spring,model_p1_spring,model_p5_spring,model_p25_spring,model_p50_spring,model_p75_spring,model_p95_spring,model_p99_spring,
model_diurnal_mag_summer,model_diurnal_ph_summer,model_mean_summer,model_p1_summer,model_p5_summer,model_p25_summer,model_p50_summer,model_p75_summer,model_p95_summer,model_p99_summer,
model_diurnal_mag_autumn,model_diurnal_ph_autumn,model_mean_autumn,model_p1_autumn,model_p5_autumn,model_p25_autumn,model_p50_autumn,model_p75_autumn,model_p95_autumn,model_p99_autumn,
model_diurnal_mag_winter,model_diurnal_ph_winter,model_mean_winter,model_p1_winter,model_p5_winter,model_p25_winter,model_p50_winter,model_p75_winter,model_p95_winter,model_p99_winter,
model_seasonal_mag_day,model_seasonal_ph_day,model_mean_day,model_p1_day,model_p5_day,model_p25_day,model_p50_day,model_p75_day,model_p95_day,model_p99_day,
model_seasonal_mag_night,model_seasonal_ph_night,model_mean_night,model_p1_night,model_p5_night,model_p25_night,model_p50_night,model_p75_night,model_p95_night,model_p99_night,
model_daily_h3_mag,model_daily_h2_mag,model_daily_h1_mag,model_daily_mag,model_annual_h3_mag,model_annual_h2_mag,model_annual_h1_mag,model_annual_mag]

#xy_params_text
xy_params_text = ['obs_diurnal_mag','obs_diurnal_ph','obs_seasonal_mag','obs_seasonal_ph','obs_mean','obs_p1','obs_p5','obs_p25','obs_p50','obs_p75','obs_p95','obs_p99','obs_pc_var_daily','obs_pc_var_seasonal','obs_pc_var_full',
'obs_diurnal_mag_spring','obs_diurnal_ph_spring','obs_mean_spring','obs_p1_spring','obs_p5_spring','obs_p25_spring','obs_p50_spring','obs_p75_spring','obs_p95_spring','obs_p99_spring',
'obs_diurnal_mag_summer','obs_diurnal_ph_summer','obs_mean_summer','obs_p1_summer','obs_p5_summer','obs_p25_summer','obs_p50_summer','obs_p75_summer','obs_p95_summer','obs_p99_summer',
'obs_diurnal_mag_autumn','obs_diurnal_ph_autumn','obs_mean_autumn','obs_p1_autumn','obs_p5_autumn','obs_p25_autumn','obs_p50_autumn','obs_p75_autumn','obs_p95_autumn','obs_p99_autumn',
'obs_diurnal_mag_winter','obs_diurnal_ph_winter','obs_mean_winter','obs_p1_winter','obs_p5_winter','obs_p25_winter','obs_p50_winter','obs_p75_winter','obs_p95_winter','obs_p99_winter',
'obs_seasonal_mag_day','obs_seasonal_ph_day','obs_mean_day','obs_p1_day','obs_p5_day','obs_p25_day','obs_p50_day','obs_p75_day','obs_p95_day','obs_p99_day',
'obs_seasonal_mag_night','obs_seasonal_ph_night','obs_mean_night','obs_p1_night','obs_p5_night','obs_p25_night','obs_p50_night','obs_p75_night','obs_p95_night','obs_p99_night',
'obs_daily_h3_mag','obs_daily_h2_mag,obs_daily_h1_mag','obs_daily_mag','obs_annual_h3_mag','obs_annual_h2_mag','obs_annual_h1_mag','obs_annual_mag',
'model_diurnal_mag','model_diurnal_ph','model_seasonal_mag','model_seasonal_ph','model_mean','model_p1','model_p5','model_p25','model_p50','model_p75','model_p95','model_p99','model_pc_var_daily','model_pc_var_seasonal','model_pc_var_full',
'model_diurnal_mag_spring','model_diurnal_ph_spring','model_mean_spring','model_p1_spring','model_p5_spring','model_p25_spring','model_p50_spring','model_p75_spring','model_p95_spring','model_p99_spring',
'model_diurnal_mag_summer','model_diurnal_ph_summer','model_mean_summer','model_p1_summer','model_p5_summer','model_p25_summer','model_p50_summer','model_p75_summer','model_p95_summer','model_p99_summer',
'model_diurnal_mag_autumn','model_diurnal_ph_autumn','model_mean_autumn','model_p1_autumn','model_p5_autumn','model_p25_autumn','model_p50_autumn','model_p75_autumn','model_p95_autumn','model_p99_autumn',
'model_diurnal_mag_winter','model_diurnal_ph_winter','model_mean_winter','model_p1_winter','model_p5_winter','model_p25_winter','model_p50_winter','model_p75_winter','model_p95_winter','model_p99_winter',
'model_seasonal_mag_day','model_seasonal_ph_day','model_mean_day','model_p1_day','model_p5_day','model_p25_day','model_p50_day','model_p75_day','model_p95_day','model_p99_day',
'model_seasonal_mag_night','model_seasonal_ph_night','model_mean_night','model_p1_night','model_p5_night','model_p25_night','model_p50_night','model_p75_night','model_p95_night','model_p99_night',
'model_daily_h3_mag','model_daily_h2_mag','model_daily_h1_mag','model_daily_mag','model_annual_h3_mag','model_annual_h2_mag','model_annual_h1_mag','model_annual_mag']


#----------------------------------------------------------------------------
#SETUP PLOT

#SETUP CALLBACKS
def onpick(event):
    ind = event.ind[0]

    #time series
    ts_i = obs_std_var[ind]
    p2.set_ydata(ts_i)
    ax2.set_ylim(np.min(ts_i)-1,np.max(ts_i)+1)
    ax2.set_title('%s Lat:%s Lon:%s Alt:%s'%(obs_refs[ind],obs_lats[ind],obs_lons[ind],obs_alt[ind]),loc='right')
    
    #diurnal waveform
    d_i = obs_diurnal_ave_waveform[ind]
    p3.set_ydata(d_i)
    ax3.set_ylim(np.min(d_i)-1,np.max(d_i)+1)

    #seasonal_waveform
    s_i = obs_seasonal_waveform[ind]
    p4.set_ydata(s_i)
    ax4.set_ylim(np.min(s_i)-1,np.max(s_i)+1)

    fig.canvas.draw()

def change_scatter_x(change):
    y_ind = xy_params_text.index(y_sel.value)
    x_ind = xy_params_text.index(change['new'])
    z_ind = xy_params_text.index(z_sel.value)
    #clear scatter plot
    ax1.cla()
    #plot new scatter
    if 'phase' not in xy_params_text[z_ind].lower():
        ax1.scatter(xy_params[x_ind],xy_params[y_ind],s=50,c=xy_params[z_ind],picker=1,vmin=np.min(xy_params[z_ind]),vmax=np.max(xy_params[z_ind]),cmap=plt.cm.gist_earth)
    else:
        z_min = 0
        if 'diurnal' in xy_params_text[z_ind].lower():
            z_max = 24
        elif 'seasonal' in xy_params_text[z_ind].lower():
            z_max = 12
        ax1.scatter(xy_params[x_ind],xy_params[y_ind],s=50,c=xy_params[z_ind],picker=1,vmin=z_min,vmax=z_max,cmap=plt.cm.hsv)
    #set axis labels
    ax1.set_xlabel(xy_params_text[x_ind])
    ax1.set_ylabel(xy_params_text[y_ind])
    #redraw
    fig.canvas.draw()

def change_scatter_y(change):
    x_ind = xy_params_text.index(x_sel.value)
    y_ind = xy_params_text.index(change['new'])
    z_ind = xy_params_text.index(z_sel.value)
    #clear scatter plot
    ax1.cla()
    #plot new scatter
    if 'phase' not in xy_params_text[z_ind].lower():
        ax1.scatter(xy_params[x_ind],xy_params[y_ind],s=50,c=xy_params[z_ind],picker=1,vmin=np.min(xy_params[z_ind]),vmax=np.max(xy_params[z_ind]),cmap=plt.cm.gist_earth)
    else:
        z_min = 0
        if 'diurnal' in xy_params_text[z_ind].lower():
            z_max = 24
        elif 'seasonal' in xy_params_text[z_ind].lower():
            z_max = 12
        ax1.scatter(xy_params[x_ind],xy_params[y_ind],s=50,c=xy_params[z_ind],picker=1,vmin=z_min,vmax=z_max,cmap=plt.cm.hsv)
    #set axis labels
    ax1.set_xlabel(xy_params_text[x_ind])
    ax1.set_ylabel(xy_params_text[y_ind])
    #redraw
    fig.canvas.draw()

def change_cb_z(change):
    x_ind = xy_params_text.index(x_sel.value)
    y_ind = xy_params_text.index(y_sel.value)
    new_z_param = change['new']
    new_z_ind = xy_params_text.index(new_z_param)
    #clear scatter plot
    ax1.cla()
    #plot new scatter
    if 'phase' not in xy_params_text[new_z_ind].lower():
        ax1.scatter(xy_params[x_ind],xy_params[y_ind],s=50,c=xy_params[new_z_ind],picker=1,vmin=np.min(xy_params[new_z_ind]),vmax=np.max(xy_params[new_z_ind]),cmap=plt.cm.gist_earth)
        #change colorbar range
        cb.set_clim(vmin=np.min(xy_params[new_z_ind]),vmax=np.max(xy_params[new_z_ind])) 
        #set colormap
        cb.set_cmap(plt.cm.gist_earth)
    else:
        z_min = 0
        if 'diurnal' in xy_params_text[new_z_ind].lower():
            z_max = 24
            #change colorbar range
            cb.set_clim(vmin=z_min,vmax=z_max)
        elif 'seasonal' in xy_params_text[new_z_ind].lower():
            z_max = 12
            #change colorbar range
            cb.set_clim(vmin=z_min,vmax=z_max)
        ax1.scatter(xy_params[x_ind],xy_params[y_ind],s=50,c=xy_params[new_z_ind],picker=1,vmin=z_min,vmax=z_max,cmap=plt.cm.hsv)
        #set colormap
        cb.set_cmap(plt.cm.hsv)
                    
    #set colorbar label
    cb.set_label(xy_params_text[new_z_ind])
    #set axis labels
    ax1.set_xlabel(xy_params_text[x_ind])
    ax1.set_ylabel(xy_params_text[y_ind])
    #redraw
    cb.draw_all() 
    fig.canvas.draw()
    
#setup dropdown boxes
#x_sel = Dropdown(description='X Axis',options=xy_params_text)
#x_sel.observe(change_scatter_x, names="value")
#y_sel = Dropdown(description='Y Axis',options=xy_params_text)
#y_sel.observe(change_scatter_y, names="value")
#z_sel = Dropdown(description='Z Axis',options=xy_params_text)
#z_sel.observe(change_cb_z, names="value")
#display(HBox((x_sel,y_sel,z_sel)))

#setup figure
fig =plt.figure(figsize=(11.5,8.5))
fig.patch.set_facecolor('white')
gs1 = gridspec.GridSpec(4, 4)
gs1.update(left=0.06,right=0.47,bottom=0.45,wspace=0.03)
ax1 = plt.subplot(gs1[:,:])
gs2 = gridspec.GridSpec(4, 4)
gs2.update(left=0.06,top=0.38,wspace=0.03)
ax2 = plt.subplot(gs2[:,:]) 
gs3 = gridspec.GridSpec(4, 4)
gs3.update(left=0.53,bottom=0.68,wspace=0.03)
ax3 = plt.subplot(gs3[:,:])         
gs4 = gridspec.GridSpec(4, 4)
gs4.update(left=0.53,top=0.65,bottom=0.45,wspace=0.03)
ax4 = plt.subplot(gs4[:,:])             
             
#xy plot
p1 = ax1.scatter(xy_params[0],xy_params[0],s=50,c=xy_params[0],picker=1,cmap=plt.cm.gist_earth,vmin=np.min(xy_params[0]),vmax=np.max(xy_params[0]))
ax1.set_xlabel(xy_params_text[0])
ax1.set_ylabel(xy_params_text[0])

#ts plot
p2, = ax2.plot(obs_datetime_time,[np.NaN]*len(time))
#ax2.set_xlim(0,50000)
             
#diurnal periodic plot
p3, = ax3.plot(obs_datetime_time[:24],[np.NaN]*len(d_time))
#ax3.set_xlim(0,24)
             
#seasonal periodic plot
p4, = ax4.plot(obs_datetime_time[:8766],[np.NaN]*len(s_time))
#ax4.set_xlim(0,8766)

#cb
cbar_ax = fig.add_axes([0.06, 0.96, 0.41, 0.04])
cb = fig.colorbar(p1, orientation = 'horizontal', cax=cbar_ax,label=xy_params_text[0])

fig.canvas.mpl_connect('pick_event', onpick)

plt.show()


 

