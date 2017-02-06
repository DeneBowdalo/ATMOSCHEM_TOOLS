from netCDF4 import Dataset
import matplotlib.pyplot as plt
import numpy as np
import os
from netCDF4 import num2date, date2num
import pandas as pd
import datetime
from netCDF4 import num2date, date2num

species = os.getcwd().split('/')[-2]

start = datetime.datetime(year = 1970, month = 1, day = 1, hour = 0, minute = 0)
end = datetime.datetime(year = 2015, month = 1, day = 1, hour = 0, minute = 0)          
grid_dates = [d.strftime('%Y%m%d') for d in pd.date_range(start,end,freq='H')][:-1]
all_years = np.array([i[:4] for i in grid_dates])

year_range = range(1970,2015)

a = Dataset('GLOBAL_SURFACE_%s_1970_2015_H_HDMN.nc'%(species))
time = a.variables['time'][:]
time = num2date(time, units='hours since 0001-01-01 00:00:00', calendar='gregorian')
refs = a.groups.keys()
for i in range(len(refs)):
    data = a.groups[refs[i]].variables[species.lower()][:]
    data[data < 0] = np.NaN
    for j in range(len(year_range)):
        year_test = all_years == str(year_range[j])
        yeardata = data[year_test]
        if np.nanmin(yeardata) > 0.09:
            print refs[i],year_range[j]                                                                                                                                    
            plt.plot(time,data)
            plt.show()
            break
