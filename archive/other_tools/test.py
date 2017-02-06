import modules
import numpy as np
from netCDF4 import Dataset
from cvxopt import solvers, matrix, spmatrix, mul
from pickle import load
solvers.options['show_progress'] = 0

#data = load(open('cvxfit.bin','rb'))
#u, y = data['u'], data['y']

#read obs netcdf file
obs_root_grp = Dataset('/work/home/db876/observations/surface/O3/process/GLOBAL_SURFACE_O3_2006_2012.nc')

#read in specific site data
site_group = obs_root_grp.groups['cvo']

#read in variables for site
obs_var = site_group.variables['o3'][:]
obs_date = site_group.variables['date'][:]
obs_time = site_group.variables['time'][:]
obs_lat = site_group.latitude
obs_lon = site_group.longitude
obs_alt = site_group.altitude
obs_group = site_group.process_group

test = obs_var >= 0 
obs_var = obs_var[test]
obs_time = modules.date_process(obs_date,obs_time,2005)
obs_time = obs_time[test]

obs_time = obs_time[:100]
obs_var = obs_var[:100]

u = np.copy(obs_time)
y = np.copy(obs_var)

m = len(u)

# minimize    (1/2) * || yhat - y ||_2^2
# subject to  yhat[j] >= yhat[i] + g[i]' * (u[j] - u[i]), j, i = 0,...,m-1
#
# Variables  yhat (m), g (m).

nvars = 2*m
P = spmatrix(1.0, range(m), range(m), (nvars, nvars))
q = matrix(0.0, (nvars,1))
q[:m] = -y

print 'stage 1'

# m blocks (i = 0,...,m-1) of linear inequalities
#
#     yhat[i] + g[i]' * (u[j] - u[i]) <= yhat[j], j = 0,...,m-1.

G = spmatrix([],[],[], (m**2, nvars))
I = spmatrix(1.0, range(m), range(m))
for i in range(m):
    print i
    # coefficients of yhat[i]
    G[list(range(i*m, (i+1)*m)), i] = 1.0

    # coefficients of g[i]
    G[list(range(i*m, (i+1)*m)), m+i] = u - u[i]

    # coefficients of yhat[j]
    G[list(range(i*m, (i+1)*m)), list(range(m))] -= I

print 'stage 2'

h = matrix(0.0, (m**2,1))

sol = solvers.qp(P, q, G, h)
yhat = sol['x'][:m]
g = sol['x'][m:]

nopts = 1000
ts = [ 2.2/nopts * t for t in range(1000) ]
f = [ max(yhat + mul(g, t-u)) for t in ts ]

try: import pylab
except ImportError: pass
else:
    pylab.figure(1, facecolor='w')
    pylab.plot(u, y, 'wo', markeredgecolor='b')
    pylab.plot(ts, f, '-g')
    pylab.axis([-0.1, 2.3, -1.1, 7.2])
    pylab.axis('off')
    pylab.title('Least-squares fit of convex function (fig. 6.24)')
    pylab.show()
