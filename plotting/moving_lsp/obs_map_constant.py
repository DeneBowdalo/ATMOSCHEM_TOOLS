#import matplotlib
#matplotlib.use('Agg')
import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib.patches import Polygon
from matplotlib.colors import LogNorm
from cmath import *
import colorsys
from scipy import stats
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker
import modules
from collections import OrderedDict
import operator
from netCDF4 import Dataset
from scipy.odr import Model, RealData, ODR, Data
from pylab import *
import pandas as pd
import matplotlib.dates as dates
from collections import Counter

fig, axes = plt.subplots(nrows=1, ncols=1,figsize=(17,10))
fig.patch.set_facecolor('white')

obs_refs = np.load('obs_refs.npy')
obs_lats = np.load('obs_lats.npy')
obs_lons = np.load('obs_lons.npy')
obs_country = np.load('obs_country.npy')
tags = modules.get_tags(obs_refs)

area_colors,area_boundaries,area_countries,area_labels = modules.area_dicts()

m = Basemap(projection='moll',lon_0 = 0,resolution='c')
                                                                                                                                                                             
m.fillcontinents(color='#cc9966',lake_color='#99ffff',alpha=0.7,zorder=0.5)
m.drawmapboundary(fill_color='#99ffff')

for i in range(len(obs_refs)):
    obs_ref = obs_refs[i]
    obs_lat = obs_lats[i]
    obs_lon = obs_lons[i]
    tag = tags[i]
    country = obs_country[i]

    area = modules.get_area(area_boundaries,country,tag,obs_lat,obs_lon,obs_ref)
  
    x, y = np.meshgrid(*m(obs_lon, obs_lat))
    m.scatter(x, y, color=area_colors[area],s=35,zorder=1,alpha=1.)

plt.tight_layout()
plt.show()
