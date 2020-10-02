import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(0, 10, 10000)
dx = x[1]-x[0]

def f(x):
    return np.sin(1*x)

def fp(x):
    return 1*np.cos(1*x)

y_FD = []
for i in range(len(x)-1):
    y_FD.append((f(x[i+1]) - f(x[i]))/dx)

y_FD.append(np.nan)
plt.plot(x,fp(x))
plt.plot(x,y_FD)
plt.show()
