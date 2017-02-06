import numpy as np
import csv
import glob
import datetime as datetime 
import matplotlib.pyplot as plt
import logging as log
import nappy
import scipy.stats as stats
import lomb_edit
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq, rfft

#reads in species name and processes unique variables
def plot(species):

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
        

    elif species == 'WINDSPEED': #Wind Speed extirpolated from UWND and VWND 
        def read_diff_species():    
            k=names.index('GMAO_UWND')
            i=names.index('GMAO_VWND')
            model_cut=np.sqrt((model[:,k]**2)+(model[:,i]**2))
            return model_cut 
        def read_diff_species_2():
            k=names2.index('GMAO_UWND')
            i=names2.index('GMAO_VWND')
            model_cut_2=np.sqrt((model2[:,k]**2)+(model2[:,i]**2))
            return model_cut_2
        units = r'$ms^{-1}$'
        first_label_pos = 1
        obs_data_name = 'Wind Speed (m/s) Campbell_(Mean)'
        unit_cut= 1
        species_type = 'Wind Speed'
        model_cut_switch = 1
        actual_species_name = 'Surface Windspeed'
       
#reads in the model data

    def readfile(filename, location):

        for files in filename:
            print files
            reader=csv.reader(open(files,'rb'), delimiter=' ', skipinitialspace = True)
            for row in reader:
                if row[1] == location: 
                    new=row[2:]
              
                    try:    
                        big.append(new)
                    
                    except:
                        big=[new]
            
                if row[0] == 'POINT':
                    names = row[2:]
   
        big=np.array(big)
        big=np.float64(big)

  
                 
        return big, names

#do I need to read everything in 
    try: 
        names
    except NameError:
# Readin the model output

        model , names = readfile(glob.glob("geos_v90103_4x5_logs/plane.log.20*"),'001') #This is oldest model version/or lower resolution
        model2, names2= readfile(glob.glob("plane.log.20*"),'CVO')# This is newer model version/ or higher resolution

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


        year=(model2[:,0]//10000)
        month=((model2[:,0]-year*10000)//100)
        day=(model2[:,0]-year*10000-month*100)

        hour=model2[:,1]//100
        min=(model2[:,1]-hour*100)

        doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(2006,1,1,0,0,0) \
              for i in range(len(year))]

        since2006_2=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]

#Define sampling intervals
    samp_spacing = 1./24.

#Convert model time array into numpy array
    since2006=np.array(since2006)
    since2006_2=np.array(since2006_2)
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

#Need to normalise model2 data also
    if model_cut_switch == 0:
        k=names2.index(species)
        model_cut_2 = model2[:,k]*unit_cut
    if model_cut_switch == 1:
        model_cut_2 = read_diff_species_2()
    standard_deviation_model_p_2 = np.std(model_cut_2)
    mean_model_p_2 = np.mean(model_cut_2)
    normal_model_2 = model_cut_2-mean_model_p_2
    normal_model_2 = normal_model_2/standard_deviation_model_p_2

#Calculate variance of pre-processed model data- should be 1 if normal
    #standard_dev_model = np.std(normal_model, dtype=np.float64)
    #variance_model = standard_dev_model**2
    #print 'Variance - pre-processed model data= ', variance_model


#Plot them all up. 
    fig=plt.figure(figsize=(20,12))

#Plot up standard conc. v time plot
    ax1= fig.add_subplot(2, 1, 1)
    fig.subplots_adjust(hspace=0.3)
    plt.plot(since2006, model_cut, color='black', label='GEOS v9.01.03 4x5 ')
    plt.plot(since2006_2, model_cut_2, color='green', label='GEOS v9.01.03 2x2.5 ')
    plt.grid(True)
    leg=plt.legend(loc=first_label_pos)
    leg.get_frame().set_alpha(0.4)
    plt.xlabel('Time (Days since 2006)')
    print units
    plt.ylabel('%s (%s)' % (species_type,units))
    plt.title('%s V Time' % (actual_species_name))

#Define sampling frequency
    samp_freq = 24

#Lomb-scargle plot
    ax3= fig.add_subplot(2, 1, 2)

#Plot axis period lines and labels
    annotate_line_y=np.arange(1e-10,1e4,1)
    freq_year = [345]*len(annotate_line_y)
    plt.plot(freq_year, annotate_line_y,'r--',alpha=0.4)
    plt.text(345, 1e-10, '1 Year', fontweight='bold')

#Obs lomb
    fa, fb, nout, jmax, prob = lomb_edit.fasper(since2006, normal_model, ofac, samp_freq)
#Divide output by sampling frequency
    fb = fb/samp_freq

#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
    #nyquist_freq_lomb_obs = frequencies[-1]
    #Si_lomb_obs = np.mean(fb)*nyquist_freq_lomb_obs
    #print nyquist_freq_lomb_obs, Si_lomb_obs, Si_lomb_obs*2 

#plot up
    plt.loglog(1./fa, fb,'kx',markersize=2, label='GEOS v9.01.03 4x5 ')

#Model lomb
    fx, fy, nout, jmax, prob2 = lomb_edit.fasper(since2006_2,normal_model_2, ofac, samp_freq)
#Divide output by sampling frequency
    fy = fy/samp_freq

#Calculate Nyquist frequency, Si and Si x 2 for normalisation checks.
    #nyquist_freq_lomb_model = frequencies[-1]
    #Si_lomb_model = np.mean(fy)*nyquist_freq_lomb_model
    #print nyquist_freq_lomb_model, Si_lomb_model, Si_lomb_model*2 

#plot up
    plt.loglog(1./fx, fy, 'gx', alpha = 0.75,markersize=2, label='GEOS v9.01.03 2x2.5 ')

    plt.grid(True)
    leg=plt.legend(loc=7)
    leg.get_frame().set_alpha(0.4)
    plt.ylim(1e-10,1e4)
    plt.xlabel('Period (Days)')
    plt.ylabel(r'PSD $(ppb^{2}/days^{-1})$')
    plt.title('Lomb-Scargle %s Power V Period' % actual_species_name)

#plt.savefig('O3_capeverde_comparison_plots.ps', dpi = 200)
    plt.show()
