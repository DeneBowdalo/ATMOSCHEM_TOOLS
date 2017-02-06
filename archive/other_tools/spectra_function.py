import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit
from lmfit import *

def spectra_f_orig(x,grad1,grad2,grad3,y_inter):
    #test = x < bp1
    #x1 = x[test]
    #test = (x >= bp1) & (x < bp2)
    #x2 = x[test]
    #test = x >= bp2
    #x3 = x[test]

    #print bp1,bp2

    x1_r = np.arange(x[0])
    x2_r = np.arange(1,x[1]+1)
    x3_r = np.arange(1,x[2]+1)

    p1 = grad1*x1_r+y_inter
    p2 = grad2*x2_r+p1[-1]
    p3 = grad3*x3_r+p2[-1]

    p = [p1]+[p2]+[p3]
    p = [item for sublist in p for item in sublist]

    return p

x = np.arange(200)
grad1 = 4
grad2 = 2
grad3 = 0.5
y_inter = 4

x_b = [70,70,60]

y = spectra_f_orig(x_b,grad1,grad2,grad3,y_inter)

#y_alt = np.copy(y)
y = y + np.random.normal(0,80,200)

#plt.plot(x,y)
#plt.show()


# fit_params = Parameters()
# fit_params['grad1'] = Parameter(value=2,min=-2,max=10)
# fit_params['grad2'] = Parameter(value=1,min=-10,max=10)
# fit_params['grad3'] = Parameter(value = 0.5,min=-10,max=10)
# fit_params['y_intercept'] = Parameter(value = 0,min=-10,max=10)
# fit_params['bp1'] = Parameter(value = 300,min=100,max=500)
# fit_params['bp2'] = Parameter(value = 800,min=600,max=1000)
# 
# result = minimize(spectra_f, fit_params, args=(y,))

#fit = spectra_f(fit_params, x,y)

#print fit

#p0 = [4,2,0.5,4,400,800]


#n_grad1 = popt[0]
#n_grad2 = popt[1]
#n_grad3 = popt[2]
#n_yinter = popt[3]

breakpoint_test = (x >= 30.0) & (x < 90.0)
bp_periods_1 = x[breakpoint_test]
bp1 = min(range(len(x)), key=lambda i: abs(x[i]-bp_periods_1[0]))

breakpoint_test = (x >= 100.0) & (x < 160.0)
bp_periods_2 = x[breakpoint_test]
bp2 = min(range(len(x)), key=lambda i: abs(x[i]-bp_periods_2[0]))

all_error = []
all_bp1 = []
all_bp2 = []
all_grad1 = []
all_grad2 = []
all_grad3 = []
all_yinter = []

for i in range(len(bp_periods_2)):
    for j in range(len(bp_periods_1)):
        bp1 = bp_periods_1[j]
        bp2 = bp_periods_2[i]
        
        test = x < bp1
        len1 = len(x[test])
        test = (x >= bp1) & (x < bp2)
        len2 = len(x[test])
        test = x >= bp2
        len3 = len(x[test])
        
        
        print bp1,bp2,len1,len2,len3
        
        x_b = [len1,len2,len3]
        
        popt, pcov = curve_fit(spectra_f_orig,x_b,y)
        print popt,pcov
        
        n_grad1 = popt[0]
        n_grad2 = popt[1]
        n_grad3 = popt[2]
        n_yinter = popt[3]
        
        y_est = spectra_f_orig(x_b,n_grad1,n_grad2,n_grad3,n_yinter)
        residuals = y - y_est
        fres = sum(residuals**2)
        all_error.append(fres)
        all_bp1.append(bp1)
        all_bp2.append(bp2)
        all_grad1.append(n_grad1)
        all_grad2.append(n_grad2)
        all_grad3.append(n_grad3)
        all_yinter.append(n_yinter)
        
bp1 = all_bp1[np.argmin(all_error)]
bp2 = all_bp2[np.argmin(all_error)]
n_grad1 = all_grad1[np.argmin(all_error)]
n_grad2 = all_grad2[np.argmin(all_error)]
n_grad3 = all_grad3[np.argmin(all_error)]
n_yinter = all_yinter[np.argmin(all_error)]

print bp1,bp2,n_grad1,n_grad2,n_grad3,n_yinter

test = x < bp1
len1 = len(x[test])
test = (x >= bp1) & (x < bp2)
len2 = len(x[test])
test = x >= bp2
len3 = len(x[test])

x_b = [len1,len2,len3]

y_est = spectra_f_orig(x_b,n_grad1,n_grad2,n_grad3,n_yinter)

plt.plot(x,y)
plt.plot(x,y_est)
plt.show()

x_b = [70,70,60]
popt, pcov = curve_fit(spectra_f_orig,x_b,y)

n_grad1 = popt[0]
n_grad2 = popt[1]
n_grad3 = popt[2]
n_yinter = popt[3]

y_est = spectra_f_orig(x_b,n_grad1,n_grad2,n_grad3,n_yinter)
residuals = y - y_est
fres = sum(residuals**2)
print fres 
