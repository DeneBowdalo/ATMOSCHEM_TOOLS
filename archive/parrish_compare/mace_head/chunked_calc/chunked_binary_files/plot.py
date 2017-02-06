import numpy as np
import glob
import string
import os
import sys
import matplotlib.pyplot as plt
from scipy.fftpack import fft,fftfreq, rfft
#from mpl_toolkits.basemap import Basemap
#from matplotlib.patches import Polygon
from cmath import *
from math import *

fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')

array = np.load('MH_O3_SFC_2008_2013.npy')
time = array[0]
vals = array[1]

print time
print vals
plt.plot(time, vals, color = 'black', marker = 'x', alpha = 0.75,markersize=2)   
plt.xlim(7305,7670)
plt.show()
