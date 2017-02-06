from scipy.odr import Model, RealData, ODR
import random
import numpy as np
import matplotlib.pyplot as plt
from scipy import stats

def f(B, x):
    '''Linear function y = m*x + b'''
    # B is a vector of the parameters.
    # x is an array of the current x values.
    # x is in the same format as the x passed to Data or RealData.
    #
    # Return an array in the same format as y passed to Data or RealData.
    return B[0]*x + B[1]
    
x = np.arange(1000.,3000.,1.)

y = np.arange(1000.,3000.,1.)

m = 20

for i in range(len(x)):
    x[i] = x[i] +random.normalvariate(0,m)

for i in range(len(y)):
    y[i] = y[i] +random.normalvariate(0,m)

def normalize(a):
    mean = np.average(a)
    squared = (a - mean)**2
    dev = (np.sqrt(np.sum(squared)))/len(squared)
    return (a - mean) / dev 

def normalize(a):
    a = np.array(a)
    return a/np.linalg.norm(a)

#x=(x-np.min(x))/np.ptp(x)
#y=(x-np.min(y))/np.ptp(y)

#x = normalize(x)
#y = normalize(y)




linear = Model(f)

mydata = RealData(x, y)

myodr = ODR(mydata, linear, beta0=[1., 0.])

myoutput = myodr.run()

myoutput.pprint()

#print myoutput.delta
#print myoutput.eps
print myoutput.beta
print myoutput.sum_square
#print myoutput.sum_square_delta
#print myoutput.sum_square_eps
x_d = myoutput.xplus
y_d = myoutput.y

print myoutput.res_var/len(x)

#print y

plt.plot(x,y, linestyle = 'None', marker = 'x')

plt.plot(x_d,y_d)

plt.show()

#check linear regression
#slope, intercept, r_value, p_value, std_err = stats.linregress(x,y)

#print slope, intercept
