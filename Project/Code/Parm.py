import numpy as np
import scipy.constants as scic
import matplotlib.pyplot as plt

R0 = 555e6;
const = (3.0/4.0*np.pi*scic.G*R0**2)**2;
print const;
def f(A,C):
    a = np.arange(0.0,100.0,10.0);
    return A*a**6-const*a**2+C;

def check(ff):
    if np.amin(ff) > 0.0:
        return True;
    else:
        #print 'no';
        return False;

parmA = [];
parmC = [];
AA = np.linspace(0,2.36e11,100);
CC = np.linspace(0,1e20,100);

for i in range(len(AA)):
    for j in range(len(CC)):
        if check(f(9307.0**2,CC[j])) == True :
            parmA.append(AA[i]);
            parmC.append(CC[j]);
print 9307.0**2;
print parmC;
print parmA[0], parmC[0];

plt.scatter(parmA,parmC);
#plt.ylim(-2e2,21e2);
#plt.xlim(48e2,6e3);
plt.ylabel('C');
plt.xlabel('A');
#plt.yscale('log');
plt.show();