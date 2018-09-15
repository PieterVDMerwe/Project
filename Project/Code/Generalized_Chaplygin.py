import scipy.integrate as scii
import scipy.constants as scic
import math as mt
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scio

class expan:
    def __init__(self,alpha,k,upperBound,A):
        self.alfa = alpha;
        self.k = k;
        self.a = np.arange(0.05,float(upperBound),0.001);
        self.upperBound = upperBound;
        self.A= A;
    #(5.6e-10)*  /555e6 (4*np.pi*scic.G/3)*
    def f(self,a):
        return 1.0/(np.sqrt(((self.A+a**(-3.0*(1.0+self.alfa)))**(1.0/(1.0+self.alfa)))*a**2-self.k));
    def ff(self,a):
        return (np.sqrt(((self.A+a**(-3.0*(1.0+self.alfa)))**(1.0/(1.0+self.alfa)))-self.k*a**(-2)));
    def calc(self):
        t = [];
        for i in range(len(self.a)):
            res, err = scii.quad(pos1.f,0.05,self.a[i]);
            t.append(res);
        return t;
    def acc(self):
        return -((self.a**(-3.0*(self.alfa+1.0))+1)**(1.0/(self.alfa+1)))*((self.a**(-3.0*(self.alfa+1.0))-2))

def func(x,a,b):
    return a*np.exp(b*x);
xdata = np.array([]);
ydata = np.array([]);
Lambda = 9307.0#1e-23;  
k = [-1.0, 0.0, 1.0];  
fig = plt.figure();     
ax = [];
l = 0;
for i in range(len(k)):
    pos1 = expan(1.0,k[i],1.0,Lambda**2);
    #pos2 = expan(3.0,1.0,50.0,19000000000000000.0);
    ax.append(fig.add_subplot(3,3,(i+1)));
    pos1calc = pos1.calc();
    if i < 1:
        xdata = np.concatenate([xdata,pos1calc]);
        ydata = np.concatenate([ydata,pos1.a]);
    ax[l].plot(pos1calc,pos1.a, label='integral');
    ax.append(fig.add_subplot(3,3,(i+4)));
    if i < 1:
        ax[l].set_ylabel('a');
    l = l+1;
    ax[l].plot(pos1.a,1.0/pos1.f(pos1.a));
    ax.append(fig.add_subplot(3,3,(i+7)));
    if i < 1:
        ax[l].set_ylabel(r"$\dot{a}$");
    l = l+1;
    ax[l].plot(pos1.a,pos1.acc());
    ax[l].set_xlabel('t\nk='+str(k[i]));
    if i < 1:
        ax[l].set_ylabel(r"$\frac{\ddot{a}}{a}$");
    l = l+1;
plt.tight_layout()
plt.savefig('CaplyginGasWithoutConstants.jpg',dpi=300)
plt.show();

##############################################################################3
thefile = open('xdata.txt', 'w');
for item in xdata:
    print>>thefile,item;
thefile.close();
#print xdata;
#print 'stuffies';
#print ydata;
thefile2 = open('ydata.txt', 'w');
for item in ydata:
    print>>thefile2,item;
thefile2.close();
popt, pcov = scio.curve_fit(func,xdata,ydata);
string = 'fit: a=%5.3f,b=%5.3f'%tuple(popt);
print string;
plt.plot(xdata,ydata, label='integral');
plt.plot(xdata,func(xdata,*popt),label='fit');
plt.legend();
plt.show();
#string = 'fit: a=' + str(popt[0])+'b=' + str(popt[1])+'c=' + str(popt[2])+'d=' + str(popt[3]) + 'e='+str(popt[4]);
#ax[l].plot(pos1calc,func(pos1.a,*popt),label=string);
#plt.show()