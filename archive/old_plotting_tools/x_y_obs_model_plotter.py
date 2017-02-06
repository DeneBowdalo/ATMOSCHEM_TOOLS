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
    #Vice versa with obs_cut
    model_cut_switch = 0
    obs_switch = 0
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
        obs_switch = 1

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

        model , names = readfile(glob.glob("plane.log.20*"),'CVO')

#now read in the observations

    myfile=nappy.openNAFile('York_merge_Cape_verde_1hr_R1.na')
    myfile.readData()

#ppy.openNAFile('York_merge_Cape_verde_1hr_R1.na')
    k_var1=myfile["VNAME"].index(obs_data_name)

# OK need to conver values from a list to a numpy array
 
    if obs_switch == 0:
        var1=np.array(myfile['V'][k_var1])
    elif obs_switch == 1:
        var1=np.array(myfile['V'][k_var1])+273.15

    valids1=var1 > 0
 
    var2=var1[valids1]



#Need to normalise model data also
    if model_cut_switch == 0:
        k=names.index(species)
        model_cut = model[:,k]*unit_cut  	
    if model_cut_switch == 1:
        model_cut = read_diff_species()



    model = []
    obs =[]

    for a, b in zip(var2, model_cut):
        if a > 0 and b > 0:
            model.append(b)
            obs.append(a)


#Plot them all up. 
    fig=plt.figure(figsize=(20,12))

    obs_max = max(obs)
    model_max=max(model)

    if obs_max > model_max:
        maximum = obs_max
    else:
        maximum = model_max
    
    obs_min = min(obs)
    model_min=min(model)

    if obs_min > model_min:
        minimum = model_min
    else:
        minimum = obs_min

    line = np.arange(minimum, maximum)
#Plot up standard conc. v time plot
    plt.plot(obs, model, 'x', markersize = 1)
    plt.plot(line, line, alpha = 0.7, color='red')
    plt.grid(True)
    plt.xlabel('Cape Verde Obs. %s (%s)' %(species_type, units))
    plt.ylabel('GEOS v9.01.03 2x2.5 %s (%s)' % (species_type,units))
    plt.title('X-Y Plot for %s' % (actual_species_name))



#plt.savefig('O3_capeverde_comparison_plots.ps', dpi = 200)
    plt.show()
