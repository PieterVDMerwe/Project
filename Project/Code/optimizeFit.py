import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scio

def f(xdata, a,b,c):
    return a*np.exp(-b*xdata)+c;

x = np.arange(0,100,1);
y = 2*f(x,3,2,4)+7;
plt.scatter(x,y);
popt, pcov = scio.curve_fit(f,x,y);
plt.plot(x,f(x,popt[0],popt[1],popt[2]),label=r"fit: $y=%5.3f e^{%5.3f}+%5.3f$"%tuple(popt));
plt.legend();
plt.show();