import matplotlib.pyplot as plt
import numpy as np

fig=plt.figure(figsize=(20,12))
fig.patch.set_facecolor('white')

a = np.load('annual_phases.npy')
print a 

years = np.arange(1988,2009)

print len(a)
print len(years)

#plt.ylim(0,24)
plt.plot(years, a)
plt.xlabel('Year')
plt.ylabel('Day')
plt.title('Mace Head Annual Phase Change')
plt.show()
