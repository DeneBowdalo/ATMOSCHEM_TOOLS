import matplotlib.pyplot as plt
import numpy as np
x = np.arange(0, 4, 0.2)
y = np.exp(-x)
el = 0.1 * np.abs(np.random.randn(len(y)))
plt.errorbar(x,y, yerr=e1, fmt='.-');
plt.show()
