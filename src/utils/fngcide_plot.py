from .functions_HRHR import fngcide
import matplotlib.pyplot as plt
import numpy as np
import matplotlib as mpl


omega = 1
theta = 9.9

x = np.linspace(0,1,100)
y = np.zeros(100)

for i in range(len(x)):
    y[i] = fngcide(omega,theta,x[i])

fig = plt.figure()
ax  = fig.add_subplot(111)
ax.plot(x,y,color='b')
ax.set_xlabel('Dose')
ax.set_ylabel('Effect')
ax.set_xlim(0,1)
ax.set_ylim(0,1)
ax.grid()
plt.show()