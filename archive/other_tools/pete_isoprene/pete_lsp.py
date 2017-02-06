import modules
import matplotlib.pyplot as plt
import lomb_phase
import numpy as np
from pandas import *
from matplotlib.ticker import FuncFormatter
import matplotlib.ticker

#set up plot
fig =plt.figure(figsize=(23,12.3))
fig.patch.set_facecolor('white')
ax = fig.add_subplot(1,1,1)

data =np.loadtxt('lanc-ptr-ms_bukit-atur_20080420.na',skiprows=35)

t = data[:,0]
x = data[:,1]

new_t = np.arange(t[0],t[-1]+0.01,0.01)

new_x = np.array([9999]*len(new_t))


inds = np.searchsorted(new_t,t)

new_x[inds] = x

valid = new_x != 9999
new_t = new_t[valid]
new_x = new_x[valid]

avgdt =new_t[1]-new_t[0]

ofac = 4

model_periods,model_mag,model_ph,model_fr,model_fi,amp_corr = modules.take_lomb(new_t,new_x,ofac,avgdt)

def form2(x, pos):
    """ This function returns a string with 3 decimal places, given the input x"""
    return '%.2f' % x   

def form5(x, pos):
    """ This function returns a string with 3 decimal places, given the input x"""
    return '%.4f' % x

ax.loglog(model_periods,model_mag,color='black',markersize=5,marker='x') 
plt.grid(True)
xformatter = FuncFormatter(form2)
yformatter = FuncFormatter(form5)
ax.xaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
ax.yaxis.set_major_formatter(matplotlib.ticker.ScalarFormatter())
plt.gca().xaxis.set_major_formatter(FuncFormatter(xformatter))
plt.gca().yaxis.set_major_formatter(FuncFormatter(yformatter))
plt.ylabel('Amplitude (ppbv)',fontsize=30)
plt.xlabel('Period (Days)',fontsize=30)
plt.title('Lomb Scargle Periodogram of Isoprene (ppbv) \nmeasured at Bukit Atur, Malaysia, on 75m tower',fontsize=30,y=1.01)
plt.xlim(0.03,100)
ax.xaxis.set_tick_params(labelsize=28)
ax.yaxis.set_tick_params(labelsize=28)

for axis in ['top','bottom','left','right']:
  ax.spines[axis].set_linewidth(2)

plt.savefig('pete_lsp.png')
plt.show()
