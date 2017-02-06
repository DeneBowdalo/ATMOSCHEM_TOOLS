import glob
import numpy as np
import csv
import datetime
from collections import OrderedDict
from netCDF4 import Dataset
from collections import Counter
import modules
from time import sleep
import pandas as pd
from scipy import stats
import os
import multiprocessing
from netCDF4 import num2date, date2num
from tzwhere import tzwhere
import pytz
import os
import psutil
import sys
from subprocess import call
tz_root = tzwhere.tzwhere(shapely=True,forceTZ=True)

process_group = 'EPA AQS'

#chunk_read_n refers to the chunk of saved AQS data that will be read 
chunk_read_n = list(sys.argv)[1]

#data range
start_year = 1970
end_year = 2014  
 
fname_species = os.getcwd().split('/')[-2]
if fname_species[:3] == 'NO2':
    no2_type = fname_species[4:]
    species = 'NO2'
else:
    no2_type = 'na'
    species = fname_species

#OUTPUT RESOLUTION CAN BE HOURLY (H), HOURLY & DAILY (HD), HOURLY, DAILY AND MONTHLY (HDM)
output_res = 'HDM'

#set run type as serial or parallel
run_type = 'parallel'

#get molecular mass and required concentration resolution of species
data_resolution,mol_mass,aqs_code,airbase_code = modules.get_spec_specs(species)
    
#SET EXIT ERROR COUNTS AT 0
inv_nometa_count = 0
inv_anyvaliddata_count = 0
inv_nokeymeta_count = 0
inv_resolution_count = 0
inv_badmeasurementmethod_count = 0

exit_counts = np.array([inv_nometa_count,inv_anyvaliddata_count,inv_nokeymeta_count,inv_resolution_count,inv_badmeasurementmethod_count])

#SET N OBS COUNTS
n_obs_all = 0
n_obs_after_nometa = 0
n_obs_after_flagsandlod = 0
n_obs_after_duplicate = 0
n_obs_after_anyvaliddata = 0
n_obs_after_nokeymeta = 0
n_obs_after_resolution = 0
n_obs_after_badmeasurementmethod = 0

n_obs_counts = np.array([n_obs_all,n_obs_after_nometa,n_obs_after_flagsandlod,n_obs_after_duplicate,n_obs_after_anyvaliddata,n_obs_after_nokeymeta,n_obs_after_resolution,n_obs_after_badmeasurementmethod])

#set exit code output 
exit_nometa_refs = []
exit_nometa_lats = []
exit_nometa_lons = []
exit_nometa_pg = []
exit_anyvaliddata_refs = []
exit_anyvaliddata_lats = []
exit_anyvaliddata_lons = []
exit_anyvaliddata_pg = []
exit_nokeymeta_refs = []
exit_nokeymeta_lats = []
exit_nokeymeta_lons = []
exit_nokeymeta_pg = []
exit_resolution_refs = []
exit_resolution_lats = []
exit_resolution_lons = []
exit_resolution_pg = []
exit_badmeasurementmethod_refs = []
exit_badmeasurementmethod_lats = []
exit_badmeasurementmethod_lons = []
exit_badmeasurementmethod_pg = []

unknown_mm = []
unknown_mm_refs = []
unknown_local_tz = []

#check if can read in memory chunk if not create memory chunks
mem_files = glob.glob('/work/home/db876/observations/surface/%s/process/AQS_split/chunk**%s.nc'%(fname_species,output_res))
print mem_files
if len(mem_files) == 0:

    files = glob.glob('/work/home/db876/observations/surface/%s/AQS/*'%(fname_species))
    files=modules.natsorted(files)

    year_array = np.arange(start_year,end_year+1)

    valid_files = []
    for i in year_array: 
        for f in files:
            if str(i) in f:
                valid_files.append(f)
            
    print valid_files

    #create ref lists
    all_site_refs = []
    test_site_refs = []
    uniq_refs = []

    #read site info
    all_meta = np.load('/work/home/db876/observations/surface/%s/process/AQS_COREMETA.npy'%(fname_species))

    meta_refs = all_meta[:,0]
    meta_lats = all_meta[:,1]
    meta_lons = all_meta[:,2]
    meta_alts = all_meta[:,3]
    meta_class = all_meta[:,4]
    meta_sitenames = all_meta[:,5]

    del_list = []
    for i in range(len(meta_refs)):
    #if meta_ref is blank then delete metadata from lists
        if (meta_refs[i] == ''):
            del_list = np.append(del_list,i)
    meta_refs = np.delete(meta_refs,del_list)
    meta_lats = np.delete(meta_lats,del_list)
    meta_lons = np.delete(meta_lons,del_list)
    meta_alts = np.delete(meta_alts,del_list)
    meta_class = np.delete(meta_class,del_list)
    meta_sitenames = np.delete(meta_sitenames,del_list)

    meta_refs = meta_refs.tolist()

    all_refs = []
    refs = []
    valid_refs  = []
    yyyymmdd = []
    hhmm = []
    all_mm = []
    sample_len = []
    all_units = []

    vals = []
    n_dup_array = []
    site_resolutions = []

    #read data
    if species != 'ISOP':
        for fi in valid_files:
            print fi
            year_refs = []
            with open(fi, 'rb') as f:
                reader = csv.reader(f)
                counter = 0
                for row in reader:
                    if counter == 0:
                        state_index = row.index("State Code")
                        county_index = row.index("County Code")
                        ref_index = row.index("Site Num")
                        val_index = row.index("Sample Measurement")
                        lod_index = row.index("MDL")
                        date_index = row.index("Date GMT")
                        time_index = row.index("Time GMT")
                        unit_index = row.index("Units of Measure")
                        method_index = row.index("Method Name")  
                        param_index = row.index("Parameter Code")
            
                    elif counter != 0:
                        if row[param_index] == str(aqs_code):
                        #get limit of detection (if do not have then set it at 99999 so data is invalid)
                            #try:
                            #    lod = float(row[lod_index])
                            #except:
                            #    lod = 0
                                #lod = 99999
                            val = float(row[val_index])
                            #check val not empty and not > 0 & if val >= lod then process line
                            state = row[state_index] #read state
                            county = row[county_index] #read county
                            ref = row[ref_index] #read ref
                            joint_ref = state+county+ref
                            year_refs.append(joint_ref)
                            all_refs.append(joint_ref)
                            dt = row[date_index]
                            time = row[time_index]
                            yyyymmdd.append(dt[:4] + dt[5:7] + dt[8:])
                            hhmm.append(time[:2] + time[3:])
                            all_mm.append(row[method_index].replace(" ", ""))
                            all_units.append(row[unit_index])
                            if (val != '') & (val >= 0):# & (val >= lod):
                                vals.append(val)
                            else:
                                vals.append(-99999)

                    counter+=1
                year_refs = set(year_refs)
                refs = np.append(refs,[i for i in year_refs]) 

    isoprene_method_dict = {'091':'flame ionisation detection','101':'gas chromatography mass spectrometry/flame ionisation detection','120':'gas chromatography flame ionisation detection','122':'gas chromatography flame ionisation detection','123':'gas chromatography flame ionisation detection',
                            '125':'gas chromatography flame ionisation detection','126':'gas chromatography flame ionisation detection','127':'gas chromatography mass spectrometry','128':'gas chromatography flame ionisation detection','129':'gas chromatography mass spectrometry','135':'gas chromatography mass spectrometry','136':'gas chromatography mass spectrometry/flame ionisation detection',
                            '142':'gas chromatography','146':'gas chromatography flame ionisation detection','147':'gas chromatography flame ionisation detection','148':'gas chromatography fourier transform infrared spectroscopy/mass spectrometry','149':'gas chromatography flame ionisation detection/mass selective detection','150':'gas chromatography mass spectrometry','152':'gas chromatography flame ionisation detection','175':'gas chromatography mass spectrometry','177':'gas chromatography flame ionisation detection/mass selective detection','210':'gas chromatography mass spectrometry'}
                        
    isoprene_samplelen_dict = {'1':1,'B':3,'3':4,'N':5,'4':6,'6':12,'7':24}

    isoprene_unit_dict = {'008':'Parts per billion','078':'Parts per billion Carbon','101':'Parts per million Carbon'}

    if species == 'ISOP':
        for fi in valid_files:
            print fi
            year_refs = []
            with open(fi, 'rb') as f:
                reader = csv.reader(f,delimiter='|')
                counter = 0
                for row in reader:
                    if counter == 0:
                        state_index = row.index("State Code")
                        county_index = row.index("County Code")
                        ref_index = row.index("Site ID")
                        sample_len_index = row.index("Sample Duration")
                        val_index = row.index("Sample Value")
                        date_index = row.index("Date")
                        time_index = row.index("Start Time")
                        unit_index = row.index("Unit")
                        method_index = row.index("Method")  
                        param_index = row.index("Parameter")
            
                    if counter > 1:
                        if row[param_index] == str(aqs_code):
                            if row[0] == 'RD': 
                                try:
                                    val = float(row[val_index])
                                except:
                                    val = -99999
                                #check val not empty and not > 0 and check val >= lod (0.01), if so process line
                                state = row[state_index] #read state
                                county = row[county_index] #read county
                                ref = row[ref_index] #read ref
                                joint_ref = state+county+ref
                                dt = row[date_index]
                                time = row[time_index]
                                sample_len = isoprene_samplelen_dict[row[sample_len_index]]
                                if sample_len == 1:
                                    yyyymmdd.append([dt[:4] + dt[4:6] + dt[6:]])
                                    hhmm.append([time[:2] + time[3:]])
                                    all_mm.append([isoprene_method_dict[row[method_index]]])
                                    all_units.append([isoprene_unit_dict[row[unit_index]]])
                                    if (val >= 0):
                                        vals.append([val])
                                    else:
                                        vals.append([-99999])
                                    n_dup_array.append([0])
                                    year_refs.append([joint_ref])
                                    all_refs.append([joint_ref])
                                    site_resolutions.append(['H'])
                                
                                #deal with sample lens > 1 hour
                                else:
                                    if output_res == 'H':
                                        continue
                                    else:
                                        yyyymmdd.append([dt[:4] + dt[4:6] + dt[6:]])
                                        hhmm.append([time[:2] + time[3:]])
                                        all_mm.append([isoprene_method_dict[row[method_index]]])
                                        all_units.append([isoprene_unit_dict[row[unit_index]]])
                                        if (val >= 0):
                                            vals.append([val])
                                        else:
                                            vals.append([-99999])
                                        year_refs.append([joint_ref])
                                        all_refs.append([joint_ref])
                                        site_resolutions.append(['D'])
                            
                                        current_year = int(yyyymmdd[-1][0][:4])
                                        current_month = int(yyyymmdd[-1][0][4:6])
                                        current_day = int(yyyymmdd[-1][0][6:])
                                        current_hh = int(hhmm[-1][0][:2])
                                        current_mm = int(hhmm[-1][0][2:])
                
                                        s = datetime.datetime(year = current_year, month = current_month, day = current_day, hour = current_hh, minute = current_mm)
                                        e = s+datetime.timedelta(hours=sample_len)
                                        day_dates = [d.strftime('%Y%m%d') for d in pd.date_range(s,e,freq='H')][1:-1]
                                        day_hours = [d.strftime('%H%M') for d in pd.date_range(s,e,freq='H')][1:-1]
                            
                                        yyyymmdd.append(day_dates)
                                        hhmm.append(day_hours)
                                        vals.append([vals[-1][0]]*len(day_dates))
                                        year_refs.append([joint_ref]*len(day_dates))
                                        all_refs.append([joint_ref]*len(day_dates))
                                        all_mm.append([isoprene_method_dict[row[method_index]]]*len(day_dates))
                                        all_units.append([isoprene_unit_dict[row[unit_index]]]*len(day_dates))
                                        site_resolutions.append(['D']*len(day_dates))

                                        #append to n duplicated array
                                        n_dup_array.append([0])
                                        n_dup_array.append([1]*len(day_dates))
            
                    counter+=1
                #flatten year_refs
                year_refs = [item for sublist in year_refs for item in sublist]
                year_refs = set(year_refs)
                refs = np.append(refs,[i for i in year_refs]) 


    refs = set(refs)
    valid_refs = sorted([i for i in refs])
    print 'all refs len = ', len(valid_refs)

    #flatten lists if species == isoprene
    if species == 'ISOP':
        yyyymmdd = [item for sublist in yyyymmdd for item in sublist]
        hhmm = [item for sublist in hhmm for item in sublist]
        all_mm = [item for sublist in all_mm for item in sublist]
        all_units = [item for sublist in all_units for item in sublist]
        vals = [item for sublist in vals for item in sublist]
        all_refs = [item for sublist in all_refs for item in sublist]
        n_dup_array = [item for sublist in n_dup_array for item in sublist]
        site_resolutions = [item for sublist in site_resolutions for item in sublist]


    all_refs = np.array(all_refs)
    yyyymmdd = np.array(yyyymmdd)
    hhmm = np.array(hhmm)
    vals = np.array(vals)
    all_mm = np.array(all_mm)
    all_units = np.array(all_units)
    if species == 'ISOP':
        n_dup_array = np.array(n_dup_array)
        site_resolutions = np.array(site_resolutions)    
    
    #split data into 10(+1) chunks for memory purposes (i.e not crashing out)
    #split data by refs evenly except for remainder at end
    chunk_n = len(valid_refs)/10
    valid_ref_chunks=[valid_refs[x:x+chunk_n] for x in xrange(0, len(valid_refs), chunk_n)]
    
    #save each ref in chunk
    for x in range(len(valid_ref_chunks)):
        chunk = valid_ref_chunks[x]
        root_grp_chunk = Dataset('AQS_split/chunk%s_%s.nc'%(x,output_res), 'w')
        solo_dim = root_grp_chunk.createDimension('solo_dim', 1)
        for y in range(len(chunk)):
            site_ref = chunk[y]
            print site_ref
            site_test = all_refs == site_ref
    
            site_yyyymmdd = yyyymmdd[site_test]
            site_hhmm = hhmm[site_test]
            site_vals = vals[site_test]
            mm_big = all_mm[site_test]
            site_units = all_units[site_test]
            if species == 'ISOP':
                n_dup_arr = n_dup_array[site_test]
                site_res = site_resolutions[site_test]
            else:
                n_dup_arr = np.zeros(len(site_vals))
                site_res = ['na']*len(site_yyyymmdd)
                
            #get site meta and set if have no meta
            if site_ref not in meta_refs:
                nometa = 'Yes'
            else:
                nometa = 'No'
            
            try:
                meta_index = meta_refs.index(site_ref)
                try:
                    lat = str(meta_lats[meta_index])
                except:
                    lat = 'na'
                try:
                    lon =  str(meta_lons[meta_index])
                except:
                    lon = 'na'
                try:
                    alt =  str(meta_alts[meta_index])
                except:
                    alt = 'na'
            except:
                pass
                
            unit = stats.mode(site_units)[0][0]
            raw_class_name = meta_class[meta_index]
            site_name = meta_sitenames[meta_index]
            
            time_dim = root_grp_chunk.createDimension('time_dim%s'%(y), len(site_yyyymmdd))
             
            
            #write out to group
            ref = root_grp_chunk.createGroup('%s'%(site_ref.lower()))

            #set variables
            date = ref.createVariable('site_yyyymmdd', str, ('time_dim%s'%(y),))
            time = ref.createVariable('site_hhmm', str, ('time_dim%s'%(y),))
            data = ref.createVariable('site_vals', 'f8', ('time_dim%s'%(y),))
            measurement_methods = ref.createVariable('mm_big', str, ('time_dim%s'%(y),))
            units = ref.createVariable('site_units', str, ('time_dim%s'%(y),))
            duplicate_array = ref.createVariable('n_dup_arr', 'f8', ('time_dim%s'%(y),))
            resolutions = ref.createVariable('site_res', str, ('time_dim%s'%(y),))
            latitude = ref.createVariable('lat', str, ('solo_dim',))
            longitude = ref.createVariable('lon', str, ('solo_dim',))
            altitude = ref.createVariable('alt', str, ('solo_dim',))
            u = ref.createVariable('unit', str, ('solo_dim',))
            classification = ref.createVariable('raw_class_name', str, ('solo_dim',))
            site_nombre = ref.createVariable('site_name', str, ('solo_dim',))
            no_m = ref.createVariable('nometa', str, ('solo_dim',))
            
            date[:] = np.array(site_yyyymmdd,dtype='object')
            time[:] = np.array(site_hhmm,dtype='object')
            data[:] = np.array(site_vals,dtype='object')
            measurement_methods[:] = np.array(mm_big,dtype='object')
            units[:] = np.array(site_units,dtype='object')
            duplicate_array[:] = np.array(n_dup_arr,dtype='object')
            resolutions[:] = np.array(site_res,dtype='object')
            latitude[:] = np.array(lat,dtype='object')
            longitude[:] = np.array(lon,dtype='object')
            altitude[:] = np.array(alt,dtype='object')
            u[:] = np.array(unit,dtype='object')
            classification[:] = np.array(raw_class_name,dtype='object')
            site_nombre[:] = np.array(site_name,dtype='object')
            no_m[:] = np.array(nometa,dtype='object')
        
        root_grp_chunk.close()
        
        call("qsub -q run AQS_split/run%s.pbs"%(x), shell=True)
    
    #delete all big data from memory
    
    del all_refs
    del yyyymmdd
    del hhmm
    del vals
    del all_mm
    del all_units
    del n_dup_array
    del site_resolutions
    del site_yyyymmdd
    del site_hhmm
    del site_vals
    del mm_big
    del site_units
    del n_dup_arr
    del site_res
    del all_meta
    del meta_lats
    del meta_lons
    del meta_alts
    del meta_class
    del meta_sitenames
    del lat
    del lon
    del alt
    del unit
    del raw_class_name
    del site_name
    del nometa
    
    #if making data chunks then make job finish 
    1+'a'
    
else:
    if int(chunk_read_n) == 0:
        call("qsub -q run AQS_split/run1.pbs", shell=True)
        call("qsub -q run AQS_split/run2.pbs", shell=True)
        call("qsub -q run AQS_split/run3.pbs", shell=True)
        call("qsub -q run AQS_split/run4.pbs", shell=True)
        call("qsub -q run AQS_split/run5.pbs", shell=True)
        call("qsub -q run AQS_split/run6.pbs", shell=True)
        call("qsub -q run AQS_split/run7.pbs", shell=True)
        call("qsub -q run AQS_split/run8.pbs", shell=True)
        call("qsub -q run AQS_split/run9.pbs", shell=True)
        call("qsub -q run AQS_split/run10.pbs", shell=True)
    

#load in a data chunk
read_chunk = Dataset('AQS_split/chunk%s_%s.nc'%(chunk_read_n,output_res))
valid_refs = [str(i) for i in read_chunk.groups]

a_site_yyyymmdd = []
a_site_hhmm = []
a_site_vals = []
a_mm_big = []
a_site_units = [] 
a_site_res = []
a_n_dup_arr = []
a_lat = []
a_lon = []
a_alt = []
a_unit = []
a_raw_class_name = [] 
a_site_name = []
a_no_meta = []

#read in site data from chunk
for site_ref in valid_refs:
    a_site_yyyymmdd.append(read_chunk.groups[site_ref].variables['site_yyyymmdd'][:])
    a_site_hhmm.append(read_chunk.groups[site_ref].variables['site_hhmm'][:])
    a_site_vals.append(read_chunk.groups[site_ref].variables['site_vals'][:])
    a_mm_big.append(read_chunk.groups[site_ref].variables['mm_big'][:])
    a_site_units.append(read_chunk.groups[site_ref].variables['site_units'][:])
    a_site_res.append(read_chunk.groups[site_ref].variables['site_res'][:])
    a_n_dup_arr.append(read_chunk.groups[site_ref].variables['n_dup_arr'][:])
    a_lat.append(read_chunk.groups[site_ref].variables['lat'][0])
    a_lon.append(read_chunk.groups[site_ref].variables['lon'][0])
    a_alt.append(read_chunk.groups[site_ref].variables['alt'][0])
    a_unit.append(read_chunk.groups[site_ref].variables['unit'][0])
    a_raw_class_name.append(read_chunk.groups[site_ref].variables['raw_class_name'][0])
    a_site_name.append(read_chunk.groups[site_ref].variables['site_name'][0])
    a_no_meta.append(read_chunk.groups[site_ref].variables['nometa'][0])

del read_chunk

#find n_hours and n_days between start and end date
d0 = datetime.date(start_year, 1, 1)
d1 = datetime.date(end_year+1, 1, 1)
delta = d1 - d0
n_months = ((end_year+1)-start_year)*12
n_days = delta.days
n_hours = n_days*24
n_points = n_hours

#create grid time array
start = datetime.datetime(year = start_year, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = end_year+1, month = 1, day = 1, hour = 0, minute = 0) 
grid_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='H')][:-1]
grid_times = [d.strftime('%H%M') for d in pd.date_range(start,end,freq='H')][:-1]

#make output dates and times
dt = pd.date_range(start,end,freq='H')[:-1].tolist()
output_res_times = date2num(dt, units='hours since 0001-01-01 00:00:00', calendar='gregorian')
obs_time_pd = pd.date_range(start = start,end = end, freq = 'H')[:-1]

#setup netcdf output
root_grp = Dataset('AQS_split/AQSchunk%s_SURFACE_%s_%s_%s_%s.nc'%(chunk_read_n,fname_species,start_year,end_year+1,output_res), 'w')
root_grp.description = 'Hourly Surface %s at AQS sites - Program written by Dene Bowdalo'%(fname_species)

# dimensions
root_grp.createDimension('flex', n_points)
root_grp.createDimension('flex4', None)

#save time out
times = root_grp.createVariable('time', 'f8', ('flex',))
times[:] = output_res_times

def site_iter_process(valid_refs,c):
    #set local counts
    inv_nometa = 0
    inv_anyvaliddata = 0
    inv_nokeymeta = 0
    inv_resolution = 0
    inv_badmeasurementmethod = 0
    n_all = 0
    n_after_nometa = 0
    n_after_flagsandlod = 0
    n_after_duplicate = 0
    n_after_anyvaliddata = 0
    n_after_nokeymeta = 0
    n_after_resolution = 0
    n_after_badmeasurementmethod = 0

    #set local unknown lists
    unknown_mm_list = []
    unknown_mm_refs_list = []
    unknown_local_tz_list = []

    data_valid = True
    
    site_ref = valid_refs[c]
    print 'ref = ',site_ref,c
    
    #read in site data from chunk
    site_yyyymmdd = a_site_yyyymmdd[c]
    site_hhmm = a_site_hhmm[c]
    site_vals = a_site_vals[c]
    mm_big = a_mm_big[c]
    site_units = a_site_units[c]
    site_res = a_site_res[c]
    n_dup_arr = a_n_dup_arr[c]
    lat = a_lat[c]
    lon = a_lon[c]
    alt = a_alt[c]
    unit = a_unit[c]
    raw_class_name = a_raw_class_name[c]
    site_name = a_site_name[c]
    no_meta = a_no_meta[c]
    country = 'United States'
    contact = 'aqsteam@epa.gov'
    
    print '1'
    
    try:
        lat = np.float32(lat)
    except:
        pass
    try:
        lon = np.float32(lon)
    except:
        pass  
    try:
        alt = np.float32(alt)
    except:
        pass
        

#process data for each site at a time
#for site_ref in valid_refs:
    #site_ref = valid_refs[c]
    #site_test = all_refs == site_ref
    #site_yyyymmdd = yyyymmdd[site_test]
    #site_hhmm = hhmm[site_test]
    #site_vals = vals[site_test]
    #mm_big = all_mm[site_test]
    #site_units = all_units[site_test]
    #if species == 'ISOP':
    #    n_dup_arr = n_dup_array[site_test]
    #    site_res = site_resolutions[site_test]
    #else:
    #    n_dup_arr = np.zeros(len(site_vals))
    
    #convert to ppb
    if (species == 'O3') or (species == 'NO') or (species == 'NO2') or (species == 'CO'):
        for i in range(len(site_vals)):
            if site_units[i] == 'Parts per million':
                site_vals[i] = site_vals[i]*1.e3
            elif site_units[i] == 'Parts per billion':
                pass
            else:
                print site_units[i]
                1+'a'
        
    # convert to ppb
    if species == 'ISOP':
        for i in range(len(site_vals)):
            #078 is Parts per billion Carbon, Isoprene has 5 Carbons
            if site_units[i] == 'Parts per billion Carbon':
                site_vals[i] = site_vals[i]/5.  
            #008 is Parts per billion
            elif site_units[i] == 'Parts per billion':
                pass
            #101 is Parts per million Carbon
            elif site_units[i] == 'Parts per million Carbon':
                site_vals[i] = (site_vals[i]/5.)*1.e3
            else:
                print site_units[i]
                1+'a'
               
    #add val to total obs count
    valid_hours_dup = np.sum(n_dup_arr)
    n_all += len(site_vals) - valid_hours_dup
    
    #get site meta
    #try:
    #    meta_index = meta_refs.index(site_ref)
    #    try:
    #        lat = np.float32(meta_lats[meta_index])
    #    except:
    #        lat = 'na'
    #    try:
    #        lon =  np.float32(meta_lons[meta_index])
    #    except:
    #        lon = 'na'
    #    try:
    #        alt =  np.float32(meta_alts[meta_index])
    #    except:
    #        alt = 'na'
    #except:
    #    pass
    
    #get local timezone
    try:
        local_tz_name = tz_root.tzNameAt(lat,lon,forceTZ=True)
        pytz_obj = pytz.timezone(local_tz_name)
        datetime_offset = pytz_obj.utcoffset(datetime.datetime(2000,1,1))
        if datetime_offset < datetime.timedelta(0):
            local_tz = -(24-int(datetime_offset.seconds/60/60))
        else:
            local_tz = int(datetime_offset.seconds/60/60)
    except:
        local_tz = 'na'
        print 'TIMEZONE NOT KNOWN, SITE IS %s'%(site_ref)
        unknown_local_tz_list.append(site_ref)
    #if species is ISOP set data_tz as local_tz
    if species == 'ISOP':
        data_tz = int(local_tz)
    else:
        data_tz = 0
    
    #test if site_ref in meta_refs, if not then exit
    #also test for ISOP if have local_tz
    
    if (no_meta == 'Yes') or (data_tz == 'na'):
        inv_nometa+=1
        print 'Site Invalid. No Metadata for ref'
        if no2_type == 'MOLYBDENUM':
            n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_obs_after_anyvaliddata,inv_nokeymeta,n_obs_after_nokeymeta,inv_resolution,n_obs_after_resolution,inv_badmeasurementmethod,n_obs_after_badmeasurementmethod = 0,0,0,0,0,0,0,0,0,0,0,0,0
        exit_c_list = np.array([inv_nometa,inv_anyvaliddata,inv_nokeymeta,inv_resolution,inv_badmeasurementmethod])
        n_c_list = np.array([n_all,n_after_nometa,n_after_flagsandlod,n_after_duplicate,n_after_anyvaliddata,n_after_nokeymeta,n_after_resolution,n_after_badmeasurementmethod])
        unknown_list = [unknown_mm_list,unknown_mm_refs_list,unknown_local_tz_list]
        meta = ['na','na','na','na','na','na','na','na','na','na','na','na']
        exit_r = 'nometa'
        return c,['na'],['na'],['na'],False,meta,exit_c_list,n_c_list,unknown_list,exit_r,np.zeros(1)
    

    valid_hours_dup = np.sum(n_dup_arr)
    n_after_nometa += len(site_vals) - valid_hours_dup
        
    #adjust dates and times if tz is not equal to 0 (only for ISOP)
    #use local tz calc to adjust times to UTC
    if species == 'ISOP':
        tz = int(data_tz)
        if tz != 0:
            for i in range(len(site_yyyymmdd)):
                #create datetime
                dt = datetime.datetime(int(site_yyyymmdd[i][:4]),int(site_yyyymmdd[i][4:6]),int(site_yyyymmdd[i][6:]),int(site_hhmm[i][:2]),int(site_hhmm[i][2:]))
                if tz > 0:
                    dt  = dt - datetime.timedelta(hours = int(tz))
                elif tz < 0:
                    dt  = dt + datetime.timedelta(hours = np.abs(int(tz)))
                site_yyyymmdd[i] = dt.strftime("%Y%m%d")
                site_hhmm[i] = dt.strftime("%H%M")
 
    #put vals into full grid
    date_con = np.array(site_yyyymmdd).astype(int)
    time_con = np.array(site_hhmm).astype(int)
    
    #remove data < 1970 and >= 2015
    test_inds = (date_con >= 19700101) & (date_con < 20150101)
    date_con = date_con[test_inds]
    time_con = time_con[test_inds]
    site_vals = site_vals[test_inds]
    mm_big = mm_big[test_inds]
    n_dup_arr = n_dup_arr[test_inds]
    
    #set st_big as 'continuous'
    st_big = ['continuous']*len(site_vals)
    
    #get obs valid
    test = site_vals >= 0
    valid_hours_dup = np.sum(n_dup_arr[test])
    n_obs_valid = len(site_vals[test]) - valid_hours_dup
    n_after_flagsandlod += n_obs_valid
    
    #create max possible grid
    full_data = np.empty(n_hours)
    full_data_after_flagsandlod = np.empty(n_hours)
    big_n_dup_array = np.zeros(n_hours)
    full_data[:] = -99999
    full_data_after_flagsandlod[:] = -99999
    
    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
    converted_time = modules.date_process(date_con,time_con,start_year)
    converted_time = np.round(converted_time,decimals=5)
    syn_grid_time = np.arange(0,n_days,1./24)
    syn_grid_time = np.round(syn_grid_time,decimals=5)
    raw_indices = np.searchsorted(syn_grid_time, converted_time, side='left')
    site_vals = np.array(site_vals)
    full_data_after_flagsandlod[raw_indices] = site_vals
    raw_st = np.copy(st_big)
    raw_mm = np.copy(mm_big)
    
    #test and remove duplicate and overlap points
    converted_time,site_vals,mm_big,st_big,n_dup_arr = modules.remove_duplicate_points(site_ref,converted_time,site_vals,mm_big,st_big,n_dup_arr,output_res)
    test = site_vals >= 0
    valid_hours_dup = np.sum(n_dup_arr[test])
    n_obs_valid = int(len(site_vals[test]) - valid_hours_dup)
    n_after_duplicate += n_obs_valid
    
    #find matching times between actual times and grid of times, return big array of indices of matched indices in grid
    indices = np.searchsorted(syn_grid_time, converted_time, side='left')
    full_data[indices] = site_vals 
    big_n_dup_array[indices] = n_dup_arr
    
    #unit = stats.mode(site_units)[0][0]
    #raw_class_name = meta_class[meta_index]
    #site_name = meta_sitenames[meta_index]
    #country = 'United States'
    #contact = 'aqsteam@epa.gov'
    all_tz = [data_tz]
    
    key_meta = [lat,lon,alt]
    
    #set site file resolution 
    if species != 'ISOP':
        file_res = 'H'
    else:
        #if all site resolutions are same continue then take first file_res
        all_same = all(x == site_res[0] for x in site_res)
        if all_same == True:
            file_res = site_res[0]
        else:
        #otherwise take highest frequency res as file_res 
            if 'M' in site_res:
                file_res = 'M'
            elif 'D' in site_res:
                file_res = 'D'
            else:
                file_res = 'H'
    
    #get sampling/instrument grids
    raw_st_grid,p_st_grid,p_st_grid_after_flagsandlod,raw_mm_grid,p_mm_grid,p_mm_grid_after_flagsandlod,unknown_mm_list,unknown_mm_refs_list = modules.sampling_and_instruments_by_processgroup(site_ref,process_group,species,raw_st,raw_mm,full_data_after_flagsandlod,full_data,raw_indices,unknown_mm_list,unknown_mm_refs_list,no2_type)

    print set(p_mm_grid)

    #do quality checks                                                                                                                                                                                                                                                                                                     
    data_valid,full_data,valid_hours_dup,p_st_grid,p_mm_grid,n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_after_anyvaliddata,inv_nokeymeta,n_after_nokeymeta,inv_resolution,n_after_resolution,inv_badmeasurementmethod,n_after_badmeasurementmethod,exit_r = modules.primary_quality_control(site_ref,species,file_res,no2_type,grid_dates,full_data,big_n_dup_array,valid_hours_dup,raw_st_grid,p_st_grid,p_st_grid_after_flagsandlod,raw_mm_grid,p_mm_grid,p_mm_grid_after_flagsandlod,data_resolution,n_obs_valid,key_meta,n_all,inv_nometa,n_after_nometa,n_after_flagsandlod,n_after_duplicate,inv_anyvaliddata,n_after_anyvaliddata,inv_nokeymeta,n_after_nokeymeta,inv_resolution,n_after_resolution,inv_badmeasurementmethod,n_after_badmeasurementmethod)
    if data_valid == False:
        exit_c_list = np.array([inv_nometa,inv_anyvaliddata,inv_nokeymeta,inv_resolution,inv_badmeasurementmethod])
        n_c_list = np.array([n_all,n_after_nometa,n_after_flagsandlod,n_after_duplicate,n_after_anyvaliddata,n_after_nokeymeta,n_after_resolution,n_after_badmeasurementmethod])
        unknown_list = [unknown_mm_list,unknown_mm_refs_list,unknown_local_tz_list]
        meta = [lat,lon,alt,'na','na','na','na','na','na','na','na','na']
        return c,['na'],['na'],['na'],False,meta,exit_c_list,n_c_list,unknown_list,exit_r,np.zeros(1)

    #set processed unit
    p_unit = 'pbbv'

    #pack meta
    meta = [lat,lon,alt,raw_class_name,file_res,unit,p_unit,data_tz,local_tz,site_name,country,contact]
    
    #if blank strings in meta then convert to 'na'
    for i in range(len(meta)):
        try:
            if meta[i].strip() == '':
                meta[i] = 'na'
        except:
            pass
    
    print set(raw_st_grid)
    print set(raw_mm_grid)
    print set(p_st_grid)
    print set(p_mm_grid)
    print meta
    
    
    exit_c_list = np.array([inv_nometa,inv_anyvaliddata,inv_nokeymeta,inv_resolution,inv_badmeasurementmethod])
    n_c_list = np.array([n_all,n_after_nometa,n_after_flagsandlod,n_after_duplicate,n_after_anyvaliddata,n_after_nokeymeta,n_after_resolution,n_after_badmeasurementmethod])
    print 'exit counts = ', exit_c_list
    print 'n obs counts = ', n_c_list

    unknown_list = [unknown_mm_list,unknown_mm_refs_list,unknown_local_tz_list]

    return c,full_data,p_st_grid,p_mm_grid,data_valid,meta,exit_c_list,n_c_list,unknown_list,'na',big_n_dup_array
                

if run_type == 'serial':
    for c in range(len(valid_refs)):
        c,full_data,p_st_grid,p_mm_grid,data_valid,meta,ec,nc,unknown,exit_code,n_dup = site_iter_process(valid_refs,c)
        #change ref name
        valid_refs[c] = 'aqs'+valid_refs[c]
        #add counts up
        exit_counts = exit_counts + ec
        n_obs_counts = n_obs_counts + nc
        #append to unknown lists
        unknown_mm=np.append(unknown_mm,unknown[0])
        unknown_mm_refs=np.append(unknown_mm_refs,unknown[1])
        unknown_local_tz=np.append(unknown_local_tz,unknown[2])
        
        #append exit refs
        if no2_type != 'MOLYBDENUM':
            if exit_code != 'na':
                if exit_code == 'nometa':
                    exit_nometa_refs.append(valid_refs[c])
                    exit_nometa_lats.append(str(meta[0]))
                    exit_nometa_lons.append(str(meta[1]))
                    exit_nometa_pg.append(process_group)
                elif exit_code == 'anyvaliddata':
                    exit_anyvaliddata_refs.append(valid_refs[c])
                    exit_anyvaliddata_lats.append(str(meta[0]))
                    exit_anyvaliddata_lons.append(str(meta[1]))
                    exit_anyvaliddata_pg.append(process_group)
                elif exit_code == 'nokeymeta':
                    exit_nokeymeta_refs.append(valid_refs[c])
                    exit_nokeymeta_lats.append(str(meta[0]))
                    exit_nokeymeta_lons.append(str(meta[1]))
                    exit_nokeymeta_pg.append(process_group)
                elif exit_code == 'resolution':
                    exit_resolution_refs.append(valid_refs[c])
                    exit_resolution_lats.append(str(meta[0]))
                    exit_resolution_lons.append(str(meta[1]))
                    exit_resolution_pg.append(process_group)
                elif exit_code == 'badmeasurementmethod':
                    exit_badmeasurementmethod_refs.append(valid_refs[c])
                    exit_badmeasurementmethod_lats.append(str(meta[0]))
                    exit_badmeasurementmethod_lons.append(str(meta[1]))
                    exit_badmeasurementmethod_pg.append(process_group)
        
        modules.write_out_data(valid_refs[c],process_group,root_grp,species,full_data,p_st_grid,p_mm_grid,data_valid,meta,n_dup)
            
elif run_type == 'parallel':
    if __name__ == '__main__':
        pool = multiprocessing.Pool(processes=16)
        results = [pool.apply_async(site_iter_process, (valid_refs,c)) for c in range(len(valid_refs))]
        big_array = [r.get() for r in results]
        pool.terminate()
    
    indices_array = []
    full_data_array = []
    p_st_grid_array = []
    p_mm_grid_array = []
    data_valid_array = []
    meta_array = []
    exit_counts_array = []
    n_obs_counts_array = []
    unknown_array = []
    exit_code_array = []
    n_dup_array = []
    
    for i in range(len(big_array)):
        cut = big_array[i]
        indices_array.append(cut[0])
        full_data_array.append(cut[1])
        p_st_grid_array.append(cut[2])
        p_mm_grid_array.append(cut[3])
        data_valid_array.append(cut[4])
        meta_array.append(cut[5])
        exit_counts_array.append(cut[6])
        n_obs_counts_array.append(cut[7])
        unknown_array.append(cut[8])
        exit_code_array.append(cut[9])
        n_dup_array.append(cut[10])

    indices_array = np.array(indices_array)
    full_data_array = np.array(full_data_array)
    p_st_grid_array = np.array(p_st_grid_array)
    p_mm_grid_array = np.array(p_mm_grid_array)
    data_valid_array = np.array(data_valid_array)
    meta_array = np.array(meta_array)
    exit_counts_array = np.array(exit_counts_array)
    n_obs_counts_array = np.array(n_obs_counts_array)
    unknown_array = np.array(unknown_array)
    exit_code_array = np.array(exit_code_array)
    n_dup_array = np.array(n_dup_array)
    
    #sort arrays by indices array for sanity
    full_data_array = full_data_array[indices_array]
    p_st_grid_array = p_st_grid_array[indices_array]
    p_mm_grid_array = p_mm_grid_array[indices_array]
    data_valid_array = data_valid_array[indices_array]
    meta_array = meta_array[indices_array]
    exit_counts_array = exit_counts_array[indices_array]
    n_obs_counts_array = n_obs_counts_array[indices_array]
    unknown_array = unknown_array[indices_array]
    exit_code_array = exit_code_array[indices_array]
    n_dup_array = n_dup_array[indices_array]
    
    for c in range(len(valid_refs)):    
        #change ref name
        valid_refs[c] = 'aqs'+valid_refs[c]
        
        #add counts up
        exit_counts = exit_counts + exit_counts_array[c]
        n_obs_counts = n_obs_counts + n_obs_counts_array[c]
        
        #append to unknown lists
        unknown_mm=np.append(unknown_mm,unknown_array[c][0])
        unknown_mm_refs=np.append(unknown_mm_refs,unknown_array[c][1])
        unknown_local_tz=np.append(unknown_local_tz,unknown_array[c][2])
        
        #append exit refs
        if no2_type != 'MOLYBDENUM':
            if exit_code_array[c] != 'na':
                if exit_code_array[c] == 'nometa':
                    exit_nometa_refs.append(valid_refs[c])
                    exit_nometa_lats.append(str(meta_array[c][0]))
                    exit_nometa_lons.append(str(meta_array[c][1]))
                    exit_nometa_pg.append(process_group)
                elif exit_code_array[c] == 'anyvaliddata':
                    exit_anyvaliddata_refs.append(valid_refs[c])
                    exit_anyvaliddata_lats.append(str(meta_array[c][0]))
                    exit_anyvaliddata_lons.append(str(meta_array[c][1]))
                    exit_anyvaliddata_pg.append(process_group)
                elif exit_code_array[c] == 'nokeymeta':
                    exit_nokeymeta_refs.append(valid_refs[c])
                    exit_nokeymeta_lats.append(str(meta_array[c][0]))
                    exit_nokeymeta_lons.append(str(meta_array[c][1]))
                    exit_nokeymeta_pg.append(process_group)
                elif exit_code_array[c] == 'resolution':
                    exit_resolution_refs.append(valid_refs[c])
                    exit_resolution_lats.append(str(meta_array[c][0]))
                    exit_resolution_lons.append(str(meta_array[c][1]))
                    exit_resolution_pg.append(process_group)
                elif exit_code_array[c] == 'badmeasurementmethod':
                    exit_badmeasurementmethod_refs.append(valid_refs[c])
                    exit_badmeasurementmethod_lats.append(str(meta_array[c][0]))
                    exit_badmeasurementmethod_lons.append(str(meta_array[c][1]))
                    exit_badmeasurementmethod_pg.append(process_group)
    
        modules.write_out_data(valid_refs[c],process_group,root_grp,species,full_data_array[c],p_st_grid_array[c],p_mm_grid_array[c],data_valid_array[c],meta_array[c],n_dup_array[c]) 

#save out processing stats
dim = root_grp.createDimension('flex2', 1)

i1 = root_grp.createVariable('invalid_nometa_count','int', ('flex2',))
i2 = root_grp.createVariable('invalid_anyvaliddata_count','int', ('flex2',))
i3 = root_grp.createVariable('invalid_nokeymeta_count','int', ('flex2',))
i4 = root_grp.createVariable('invalid_resolution_count','int', ('flex2',))
i5 = root_grp.createVariable('invalid_badmeasurementmethod_count','int', ('flex2',))

n1 = root_grp.createVariable('n_obs_all','int', ('flex2',))
n2 = root_grp.createVariable('n_obs_after_nometa','int', ('flex2',))
n3 = root_grp.createVariable('n_obs_after_flagsandlod','int', ('flex2',))
n4 = root_grp.createVariable('n_obs_after_duplicate','int', ('flex2',))
n5 = root_grp.createVariable('n_obs_after_anyvaliddata','int', ('flex2',))
n6 = root_grp.createVariable('n_obs_after_nokeymeta','int', ('flex2',))
n7 = root_grp.createVariable('n_obs_after_resolution','int', ('flex2',))
n8 = root_grp.createVariable('n_obs_after_badmeasurementmethod','int', ('flex2',))

#site number exit counts
i1[:] = exit_counts[0]
i2[:] = exit_counts[1]
i3[:] = exit_counts[2]
i4[:] = exit_counts[3]
i5[:] = exit_counts[4]

print exit_counts

#n obs counts after checks
n1[:] = n_obs_counts[0]
n2[:] = n_obs_counts[1]
n3[:] = n_obs_counts[2]
n4[:] = n_obs_counts[3]
n5[:] = n_obs_counts[4]
n6[:] = n_obs_counts[5]
n7[:] = n_obs_counts[6]
n8[:] = n_obs_counts[7]

print n_obs_counts

print 'Unknown mm = ', set(unknown_mm)
print 'Unknown mm Sites = ', unknown_mm_refs 
print 'Unknown local tz sites = ', unknown_local_tz

dim = root_grp.createDimension('flex3', None)
o1 = root_grp.createVariable('UNKNOWN_MM',str, ('flex3',))
o2 = root_grp.createVariable('UNKNOWN_MM_REFS',str, ('flex3',))
o3 = root_grp.createVariable('UNKNOWN_TZ_REFS',str, ('flex3',))

e1a = root_grp.createVariable('exit_nometa_refs',str, ('flex3',))
e1b = root_grp.createVariable('exit_nometa_lats',str, ('flex3',))
e1c = root_grp.createVariable('exit_nometa_lons',str, ('flex3',))
e1d = root_grp.createVariable('exit_nometa_pg',str, ('flex3',))
e2a = root_grp.createVariable('exit_anyvaliddata_refs',str, ('flex3',))
e2b = root_grp.createVariable('exit_anyvaliddata_lats',str, ('flex3',))
e2c = root_grp.createVariable('exit_anyvaliddata_lons',str, ('flex3',))
e2d = root_grp.createVariable('exit_anyvaliddata_pg',str, ('flex3',))
e3a = root_grp.createVariable('exit_nokeymeta_refs',str, ('flex3',))
e3b = root_grp.createVariable('exit_nokeymeta_lats',str, ('flex3',))
e3c = root_grp.createVariable('exit_nokeymeta_lons',str, ('flex3',))
e3d = root_grp.createVariable('exit_nokeymeta_pg',str, ('flex3',))
e4a = root_grp.createVariable('exit_resolution_refs',str, ('flex3',))
e4b = root_grp.createVariable('exit_resolution_lats',str, ('flex3',))
e4c = root_grp.createVariable('exit_resolution_lons',str, ('flex3',))
e4d = root_grp.createVariable('exit_resolution_pg',str, ('flex3',))
e5a = root_grp.createVariable('exit_badmeasurementmethod_refs',str, ('flex3',))
e5b = root_grp.createVariable('exit_badmeasurementmethod_lats',str, ('flex3',))
e5c = root_grp.createVariable('exit_badmeasurementmethod_lons',str, ('flex3',))
e5d = root_grp.createVariable('exit_badmeasurementmethod_pg',str, ('flex3',))

o1[:] = np.array(list(set(unknown_mm)),dtype='object')
o2[:] = np.array(unknown_mm_refs,dtype='object')
o3[:] = np.array(unknown_local_tz,dtype='object')

e1a[:] = np.array(exit_nometa_refs,dtype='object')
e1b[:] = np.array(exit_nometa_lats,dtype='object')
e1c[:] = np.array(exit_nometa_lons,dtype='object')
e1d[:] = np.array(exit_nometa_pg,dtype='object')
e2a[:] = np.array(exit_anyvaliddata_refs,dtype='object')
e2b[:] = np.array(exit_anyvaliddata_lats,dtype='object')
e2c[:] = np.array(exit_anyvaliddata_lons,dtype='object')
e2d[:] = np.array(exit_anyvaliddata_pg,dtype='object')
e3a[:] = np.array(exit_nokeymeta_refs,dtype='object')
e3b[:] = np.array(exit_nokeymeta_lats,dtype='object')
e3c[:] = np.array(exit_nokeymeta_lons,dtype='object')
e3d[:] = np.array(exit_nokeymeta_pg,dtype='object')
e4a[:] = np.array(exit_resolution_refs,dtype='object')
e4b[:] = np.array(exit_resolution_lats,dtype='object')
e4c[:] = np.array(exit_resolution_lons,dtype='object')
e4d[:] = np.array(exit_resolution_pg,dtype='object')
e5a[:] = np.array(exit_badmeasurementmethod_refs,dtype='object')
e5b[:] = np.array(exit_badmeasurementmethod_lats,dtype='object')
e5c[:] = np.array(exit_badmeasurementmethod_lons,dtype='object')
e5d[:] = np.array(exit_badmeasurementmethod_pg,dtype='object')

print 'NOMETA_EXIT_REFS = ', exit_nometa_refs
print 'ANYVALIDDATA_EXIT_REFS = ', exit_anyvaliddata_refs
print 'NOKEYMETA_EXIT_REFS = ', exit_nokeymeta_refs
print 'RESOLUTION_EXIT_REFS = ', exit_resolution_refs
print 'BADMEASUREMENTMETHOD_EXIT_REFS = ', exit_badmeasurementmethod_refs

print 'job finished'

root_grp.close()
