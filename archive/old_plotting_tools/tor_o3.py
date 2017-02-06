import numpy as np
import csv
import glob
import datetime as datetime 
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib import colors, ticker, cm
import logging as log
import nappy
import scipy.stats as stats
from scipy.interpolate import interp1d
from scipy.fftpack import fft,fftfreq
import datetime

#Loop through all runs, set switch to be on or off, to determine if you want to loop or not
flight_loop_switch = 'on'
species_loop_switch = 'off'

#Input flightname of interest if not looping
one_flight = 'TORERO_RF08'

#Input species of interest if not looping
one_species = 'TRA_47'


#Lists of species and flights
loop_flights=['TORERO_TF01','TORERO_TF02','TORERO_FF01','TORERO_RF01','TORERO_RF02','TORERO_RF03','TORERO_RF04','TORERO_RF05','TORERO_RF06','TORERO_RF07','TORERO_RF08','TORERO_RF09','TORERO_RF10','TORERO_RF11','TORERO_RF12','TORERO_RF13','TORERO_RF14','TORERO_RF15','TORERO_RF16','TORERO_RF17','TORERO_FF02']

species_list=['O3','CO','NO','NO2','NO3','N2O5', 'HNO4', 'HNO3', 'HNO2','PAN','PPN','PMN','R4N2','H2O2','MP','CH2O','HO2','OH','RO2','MO2','ETO2','CO','C2H6','C3H8','PRPE','ALK4','ACET','ALD2','MEK','RCHO','MVK','SO2','DMS','MSA','SO4','TRA_1','TRA_2','TRA_14','TRA_28','TRA_30','TRA_31','TRA_32','TRA_33','TRA_34','TRA_35','TRA_36','TRA_37','TRA_38','TRA_39','TRA_40','TRA_41','TRA_42','TRA_43','TRA_44','TRA_45','TRA_46','TRA_47','TRA_48','TRA_49','TRA_50','TRA_51','TRA_52','TRA_53','GMAO_TEMP','GMAO_ABSH','GMAO_PSFC','GMAO_UWND','GMAO_VWND']

tracer_list=['TRA_1','TRA_2','TRA_14','TRA_28','TRA_30','TRA_31','TRA_32','TRA_33','TRA_34','TRA_35','TRA_36','TRA_37','TRA_38','TRA_39','TRA_40','TRA_41','TRA_42','TRA_43','TRA_44','TRA_45','TRA_46','TRA_47','TRA_48','TRA_49','TRA_50','TRA_51','TRA_52','TRA_53']

tracer_names=['NOx','Ox','MACR','SO4s','NH3','NH4','NIT','NITs','BCPI','OCPI','BCPO','OCPO','DST1','DST2','DST3','DST4','SALA','SALC','Br2','Br','BrO','HOBr','HBr','BrNO2','BrNO3','CHBr3','CH2Br2','CH3Br']

ppbv_list=['O3','CO','TRA_1','TRA_2','TRA_45','TRA_46','TRA_47','TRA_48','TRA_49','TRA_50','TRA_51','TRA_52','TRA_53']

pptv_list = ['NO','NO2','NO3','N2O5','HNO4','HNO3','HNO2',' PAN','PPN','PMN','R4N2','H2O2','MP','CH2O','HO2','OH','RO2','MO2','ETO2','C2H6','C3H8','PRPE','ALK4','ACET','ALD2','MEK','RCHO','MVK','SO2','DMS','MSA','SO4','TRA_14','TRA_28','TRA_30','TRA_31','TRA_32','TRA_33','TRA_34','TRA_35','TRA_36','TRA_37','TRA_38','TRA_39','TRA_40','TRA_41','TRA_42','TRA_43','TRA_44']

def readfile(filename, location):

    read = np.load(filename)
    names = read[0,2:]
    locs = np.where(read[:,1] == location)
    big = read[locs]
    valid_data = big[:,2:]
    names = names.tolist()   
    valid_data = np.float64(valid_data)   
 
    return valid_data, names


def run(flight,species,code):
    #Assign indices to variables
    #Torero Flights

    print species
#Flight Number 1
    if flight == 'TORERO_TF01':
        flightname = 'TORERO_TF01'
        startdate  = 20120106
        enddate    = 20120106

#Flight Number 2
    if flight == 'TORERO_TF02':
        flightname = 'TORERO_TF02'
        startdate  = 20120114
        enddate    = 20120114

#Flight Number 3
    if flight == 'TORERO_FF01':
        flightname = 'TORERO_FF01'
        startdate  = 20120118
        enddate    = 20120118

#Flight Number 4
    if flight == 'TORERO_RF01':
        flightname = 'TORERO_RF01'
        startdate  = 20120119
        enddate    = 20120119

#Flight Number 5
    if flight == 'TORERO_RF02':
        flightname = 'TORERO_RF02'
        startdate  = 20120121
        enddate    = 20120121

#Flight Number 6
    if flight == 'TORERO_RF03':
        flightname = 'TORERO_RF03'
        startdate  = 20120124
        enddate    = 20120124

#Flight Number 7
    if flight == 'TORERO_RF04':
        flightname = 'TORERO_RF04'
        startdate  = 20120127
        enddate    = 20120127    

#Flight Number 8
    if flight == 'TORERO_RF05':
        flightname = 'TORERO_RF05'
        startdate  = 20120129
        enddate    = 20120129

#Flight Number 9
    if flight == 'TORERO_RF06':
        flightname = 'TORERO_RF06'
        startdate  = 20120131
        enddate    = 20120131

#Flight Number 10
    if flight == 'TORERO_RF07':
        flightname = 'TORERO_RF07'
        startdate  = 20120203
        enddate    = 20120203

#Flight Number 11
    if flight == 'TORERO_RF08':
        flightname = 'TORERO_RF08'
        startdate  = 20120204
        enddate    = 20120204

#Flight Number 12
    if flight == 'TORERO_RF09':
        flightname = 'TORERO_RF09'
        startdate  = 20120207
        enddate    = 20120207

#Flight Number 13
    if flight == 'TORERO_RF10':
        flightname = 'TORERO_RF10'
        startdate  = 20120210
        enddate    = 20120210

#Flight Number 14
    if flight == 'TORERO_RF11':
        flightname = 'TORERO_RF11'
        startdate  = 20120212
        enddate    = 20120212

#Flight Number 15
    if flight == 'TORERO_RF12':
        flightname = 'TORERO_RF12'
        startdate  = 20120214
        enddate    = 20120214

#Flight Number 16
    if flight == 'TORERO_RF13':
        flightname = 'TORERO_RF13'
        startdate  = 20120217
        enddate    = 20120217

#Flight Number 17
    if flight == 'TORERO_RF14':
        flightname = 'TORERO_RF14'
        startdate  = 20120219
        enddate    = 20120219

#Flight Number 18
    if flight == 'TORERO_RF15':
        flightname = 'TORERO_RF15'
        startdate  = 20120222
        enddate    = 20120222

#Flight Number 19
    if flight == 'TORERO_RF16':
        flightname = 'TORERO_RF16'
        startdate  = 20120224
        enddate    = 20120224

#Flight Number 20
    if flight == 'TORERO_RF17':
        flightname = 'TORERO_RF17'
        startdate  = 20120226
        enddate    = 20120226

#Flight Number 21
    if flight == 'TORERO_FF02':
        flightname = 'TORERO_FF02'
        startdate  = 20120229
        enddate    = 20120229

    if species in ppbv_list: 
        units = 'ppbV'
        var_type = 'Conc.'
        unit_convert = 1e9

    if species in pptv_list: 
        units = 'pptV'
        var_type = 'Conc.'
        unit_convert = 1e12

    if species == 'GMAO_TEMP':
        units = 'K'
        var_type = 'Temp.'
        unit_convert = 1

    if species == 'GMAO_ABSH':
        units = '%'
        var_type = 'Abs. Humid.'
        unit_convert = 1

    if species == 'GMAO_PSFC':
        units = 'hPa'
        var_type = 'Pressure'
        unit_convert = 1

    if species == 'GMAO_UWND':
        units = 'm/s'
        var_type = 'Windspeed'
        unit_convert = 1

    if species == 'GMAO_VWND':
        units = 'm/s'
        var_type = 'Windspeed'
        unit_convert = 1


    lat_index = names.index('LAT')   
    lon_index = names.index('LON')
    pre_index = tor_names.index('PRESS')
    k=names.index(species)

    if species[0:2] == 'TR':
        tracer_index =tracer_list.index(species) 
        species = tracer_names[tracer_index]      
  
    #Process valid dates to cut out a flight
    start_date_testing = model[:,0] >= startdate
    cut_model = model[start_date_testing]
    end_date_testing = cut_model[:,0] <= enddate
    cut_model = cut_model[end_date_testing]
    variable = cut_model[:,k]*unit_convert
    pressures = cut_model[:,pre_index]


    start_time = cut_model[0,1]
    end_time   = cut_model[-1,1]
   
    #Sort variable array into 2d array for contour plot
    for i in range(21):
        new_list = variable[i::21]
        new_list = np.array(new_list)
        try:
            sorted_variable =np.vstack((sorted_variable,new_list))
        except:
            sorted_variable = [new_list]
            sorted_variable=np.array(sorted_variable)
 
    
    #Get Torero flight points
    start_date_testing = tor_model[:,0] >= startdate
    cut_model_torero = tor_model[start_date_testing]
    end_date_testing = cut_model_torero[:,0] <= enddate
    cut_model_torero = cut_model_torero[end_date_testing]
    start_time_testing= cut_model_torero[:,1] >= start_time
    cut_model_torero = cut_model_torero[start_time_testing]
    end_time_testing= cut_model_torero[:,1] <= end_time
    cut_model_torero = cut_model_torero[end_time_testing]
    torero_pressure_testing = cut_model_torero[:,pre_index] > 0
    cut_model_torero = cut_model_torero[torero_pressure_testing]
    torero_pressures = cut_model_torero[:,pre_index]
     
    #Process time
    year=(cut_model[:,0]//10000)
    month=((cut_model[:,0]-year*10000)//100)
    day=(cut_model[:,0]-year*10000-month*100)
    
    hour=cut_model[:,1]//100
    min=(cut_model[:,1]-hour*100)
    
    doy=[ datetime.datetime(np.int(year[i]),np.int(month[i]),np.int(day[i]),\
                                np.int(hour[i]),np.int(min[i]),0)- \
              datetime.datetime(2012,1,1,0,0,0) \
              for i in range(len(year))]
    
    since2012=[doy[i].days+doy[i].seconds/(24.*60.*60.) for i in range(len(doy))]

    since2012=np.array(since2012)

#Crop duplicate times in array
    newsince2012, indices = np.unique(since2012,return_index=True)
    since2012 = since2012[indices]


    start_ref = since2012[0]
    end_ref = since2012[-1]

#Plot them all up. 
    fig=plt.figure(figsize=(11,9))

#Create contour grid to plot on
    xlist = since2012
    ylist = np.arange(50, 1050, 50)
    ylist = np.append(ylist, 1030)

    z = sorted_variable
    
    var_max = np.max(variable)+0.001
    var_min_test = variable > 0
    var_altered = variable[var_min_test]
    var_min = np.min(var_altered)
    
#Control color bar range
    crange = np.arange(var_min,var_max,0.001)
    plt.contourf(xlist, ylist, z, crange)
    print len(since2012)
    print len(torero_pressures)      
    plt.plot(since2012, torero_pressures, 'k:', label='Flight track')
    plt.ylim(50,1030)
    plt.gca().invert_yaxis()
    plt.colorbar()
    plt.xlabel('Time (Days since 2012)')
    plt.ylabel('Pressure, hPa' ) 
    plt.text(start_ref,70,'%s      Date:%8i' %(flightname,startdate),  fontweight='bold')
    plt.text(start_ref,90,'Start time:%i   End time:%i' %(start_time, end_time),  fontweight='bold')
    plt.text(end_ref+0.005,33,'%s %s (%s)' %(species,var_type,units))
    leg=plt.legend(loc=1)
    leg.get_frame().set_alpha(0)
    plt.title('Curtain Plot of %s %s V Time with Pressure' %(species, var_type))

    if code == 1: 
        fig.savefig('allspecies_allflights/%s_%s.pdf' %(species,flightname))
    if code == 2:
        fig.savefig('onespecies_allflights/%s_%s.pdf' %(species,flightname))
    if code == 3:
        fig.savefig('allspecies_oneflight/%s_%s.pdf' %(species,flightname))
    if code == 4:
        fig.savefig('onespecies_oneflight/%s_%s.pdf' %(species,flightname)) 
   #plt.show()



try:
    names
except NameError:
# Read in the model output

    model, names = readfile("torero_logs.npy",'COL')
    tor_model, tor_names = readfile("torero_logs.npy",'TOR')


#Loops through all flights for all species
if flight_loop_switch == 'on' and species_loop_switch == 'on':
    for i in species_list:
        for j in loop_flights:
            flight = j
            species = i
            code = 1
            run(flight,species, code)

#Loops through all flights for a certain species
if flight_loop_switch == 'on' and species_loop_switch == 'off':
    for i in loop_flights:
        flight = i
        species = one_species
        code = 2
        run(flight,species, code)
          
#Loops through all species for a certain flight
if flight_loop_switch == 'off' and species_loop_switch == 'on':
    for i in species_list:
        flight = one_flight
        species = i
        code = 3
        run(flight,species,code)

#Certain species for a certain flight   
if flight_loop_switch == 'off' and species_loop_switch == 'off':
    flight = one_flight
    species = one_species
    code = 4
    run(flight,species,code) 
