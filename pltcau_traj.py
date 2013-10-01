import numpy as np
import matplotlib.pyplot as plt

def readData(filename, col1, col2):
	data = np.loadtxt(filename)
	return data[:,col1],data[:,col2]

xcau,ycau = readData('caustics.dat',0,1)
xs,ys = readData('fort.79',1,2)

t0,mag0 = readData('data/curve-4-4225.dat',0,1)
t1,mag1 = readData('fort.35',0,1)
residual = np.mean(mag1-mag0)
mag1 = [mag1i-residual for mag1i in mag1]

plt.figure()
#plt.axis('equal')
#plt.scatter(xcau,ycau,color='r',s=1)
#plt.scatter(xs,ys,color='g',s=5)

plt.scatter(t0,mag0, color='r', s=1)
plt.scatter(t1,mag1, color='g', s=1)
#plt.scatter(t0,residual,s=1)
plt.gca().invert_yaxis()
plt.show()
