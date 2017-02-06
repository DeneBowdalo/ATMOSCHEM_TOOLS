import matplotlib.pyplot as plt
x = [5,3,7,2,4,1,11,25,33]
plt.plot(x)
plt.xticks(range(len(x)), ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i']);
plt.yticks(range(1,36,2));
plt.show()
