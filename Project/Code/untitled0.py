import scipy.integrate as scii
import scipy.constants as scic
import math as mt
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scio
from scipy.special import hyp2f1

def f(x):
    return 1.0*x - 1.0/x;

x = np.arange(0.0,10.0,0.001);
plt.plot(x,f(x));
plt.ylim(-30.0,30.0);
plt.show();