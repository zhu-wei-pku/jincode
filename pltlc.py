import numpy as np
import matplotlib.pyplot as plt

data = np.loadtxt('test-lc')
time = data[:,0]
mu = data[:,1]

plt.scatter(time,mu)
#plt.gca().invert_yaxis()
plt.show()
