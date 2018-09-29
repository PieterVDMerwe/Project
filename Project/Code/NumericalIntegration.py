import scipy.integrate as scii
import scipy.constants as scic
import math as mt
import numpy as np
import matplotlib.pyplot as plt
import scipy.optimize as scio
from scipy.special import hyp2f1

############################Integration Function###############################
def intag(func,upper ,lower):
    a = np.arange(lower,upper,0.001);
    t = [];
    for i in a:
        res, err = scii.quad(func,lower,i);
        t.append(res);
    plt.plot(t,a,label='numerical, k='+str(k));
    plt.xlabel('time t');
    plt.ylabel('scale factor a');
#    plt.show();

###############################################################################
#                   Chaplygin gas
###############################################################################
print '###############################################################################\n                 Chaplygin gas                  \n###############################################################################'; 
######################Numerical integration####################################
A=((8.0*np.pi*scic.G)/(3.0*(scic.c**2.0)))
#A1 = 1.0;
#A2 = 1.0;
#c2 = 1.0;
A1 = 1.0;
A2 = 1.0;
c2 = 1.0;
B1=1.0+A2;
B2=A1/B1;
B3=c2/B1;
beta=2.0;
k=0.0#-0.005;
F=1.0/(55.5e7)**2;

def adotParm(A,B1,B2,B3,beta,k,F,a):
    return 1.0/np.sqrt(A*((B3*a**(-3.0*beta*B1) + B2)**(1.0/beta))*(a**2) - k*F);

def adotNoParm(a):
    return adotParm(A,B1,B2,B3,beta,k,F,a);

def chz(z):
    return -((1.0+z)**(-2))*((((B3*(1.0+z)**(beta*B1)+B2*(z+1.0)**(-2.0*beta))**(1.0/beta))-k)**(-1.0/2.0));

def chadot(a):
    return np.sqrt(A*(B3*a**(-3.0*beta*B1) + B2)**(1.0/beta) - k*F*(a**2));

    
#########################Analytical solution###################################
def analyt(upper, lower):
    a = np.arange(lower,upper,0.001);
    firstFact = (1.0/((3.0/2.0)*(A**0.5)*(B2**(1.0/(2.0*beta)))*(B1)));
    arg = (B3/B2)*a**(-3.0*B1*beta) +1.0;
    t = firstFact*((arg**(-1.0/(2.0*beta)) + (1.0/(2.0*beta+1))*(arg**(-1.0 - 1.0/(2.0*beta)))*(hyp2f1(1.0,1.0+1.0/(2.0*beta),2.0+1.0/(2.0*beta),arg**(-1.0)))));
    plt.plot(t,a,label='analytical, k=0',linewidth=3.0,dashes=[2, 2]);
    

def analytChz(upper, lower):
    z = np.arange(lower,upper,0.001);
    a = 1.0/(1.0+z);
    firstFact = (1.0/((3.0/2.0)*(A**0.5)*(B2**(1.0/(2.0*beta)))*(B1)));
    arg = (B3/B2)*a**(-3.0*B1*beta) +1.0;
    t = firstFact*((arg**(-1.0/(2.0*beta)) + (1.0/(2.0*beta+1))*(arg**(-1.0 - 1.0/(2.0*beta)))*(hyp2f1(1.0,1.0+1.0/(2.0*beta),2.0+1.0/(2.0*beta),arg**(-1.0)))));
    plt.plot(t,a,label='analytical, k=0',linewidth=3.0,dashes=[2, 2]);   

def addot(a):
    rho = (B3*a**(-3.0*beta*B1) + B2)**(1.0/beta);
    P = A2*rho - A1/(rho**(beta-1.0));
    return -(rho+3.0*P)/6.0;
    
def chrho(a):
    return (B3*a**(-3.0*beta*B1) + B2)**(1.0/beta);

upper = 2.0
ak = [-1.0,0.0,1.0];
for i in ak:
    k=i;
    intag(adotNoParm,upper,0.001);
    if (k==0.0):
        analyt(upper,0.0001);
    plt.legend();
    plt.ylabel('a');
    plt.xlabel('t');
    plt.xscale('log');
plt.legend();
plt.savefig('a_ch.jpg',dpi=300);
plt.show();

upper = 2.0
ak = [-1.0,0.0,1.0];
for i in ak:
    k=i;
    intag(chz,upper,0.001);
#    if (k==0.0):
#        analytChz(upper,0.0001);
    plt.legend();
    plt.ylabel('z');
    plt.xlabel('t');
#    plt.xscale('symlog');
plt.legend();
plt.savefig('z_ch.jpg',dpi=300);
plt.show();

for i in ak:
    k=i;
    x = np. arange(0.001,10.0,0.001);
    plt.plot(x,chadot(x),label='k='+str(k));
    plt.yscale('log');
    plt.ylabel(r"$H$");
    plt.xlabel(r"a");
plt.legend();
plt.savefig('ch_H.jpg',dpi=300);
plt.show();

upper = 2.0
lw = [7.0,3.0,1.0];
for i in range(3):
    k=ak[i];
    plt.subplot(1,3,i+1);
    a = np. arange(0.0001,upper,0.00001);
    q = -A*addot(a)*(a**2)/(1.0/adotNoParm(a))**(2.0);
    plt.plot(a,q,label='k='+str(k));
#    plt.plot(a,q,label='k='+str(k),linewidth=lw[i]);
    #plt.plot(a,addot(a));
    if i == 0:
        plt.ylabel(r"$q$");
#        plt.yscale('symlog');
        plt.ylim(-0.5,1.0)
    if i == 2:
        plt.yscale('log');
    if i == 1:
        plt.xlabel('a');
    plt.grid(True);
#    plt.ylim(-2.9e18,2.0e3);
#    plt.ylim(-0.25e26,10.0e25);
    #plt.yscale('log');
    plt.legend();
plt.tight_layout();
plt.savefig('ch_q.jpg',dpi=300);
plt.show();

for i in range(3):
    k=ak[i];
    plt.plot(a,addot(a),label='k='+str(k),linewidth=lw[i]);
    #plt.plot(a,addot(a));
    plt.ylabel(r"$\ddot{a}/a$");
    plt.xlabel('a');
    plt.ylim(-2.0,0.5);
    plt.grid(True);
    #plt.yscale('log');
plt.legend();
plt.savefig('ch_ddota.jpg',dpi=300);
plt.show();

x = np. arange(0.001,15,0.1);
plt.plot(x,chrho(x));
plt.yscale('symlog');
plt.ylabel(r"$\rho$");
plt.xlabel(r"a");
plt.savefig('Ch_rho.jpg',dpi=300);
plt.show();


###############################################################################
###################         UDF             ###################################  
###############################################################################

print '###############################################################################\n                 UDF                  \n###############################################################################'; 

c1=-0.0000000000005;
c2=2.0*c1;
c3=-0.7;
c4=1.0-c1-c3;
k=-1.0;
def UDFF(a):
    return 1.0/np.sqrt((c1*a**3-c2*a-c3*a**2+c4/a)-k);

def UDFH(a):
    answ = (c1*a-c2/a-c3+c4/a**3)-k/a**2;
    #answ[answ < 0] = 0.0;  
    return answ;#np.sqrt(answ);

def UDFrho(a,Pa,Pb,C):
    return -Pa+3.0*Pb*(a-2.0/a)/4.0 +C*a**(-3.0);

def UDFacc(a,Pa,Pb,C):
    P = Pa + Pb*(a**(-1.0)-a);
    return -(UDFrho(a,Pa,Pb,C)+3.0*P)/6.0;

def UDFz(z):
    return -((1+z)**(-2))/(np.sqrt((c1*((1+z)**(-3))-c2*((1+z)**(-1))-c3*((1+z)**(-2))+c4*(1+z))-k));
    
#    return 1.0/np.sqrt(A*(c1*a**3-c2*a-c3*a**2+c4/a)-k*F);
#def fff(x):
#    return 1.0/np.sqrt(2.0*(1.0-2.0*np.exp(-2.0*x))*(2.0*np.exp(x)-1.0));
#intag(UDFF,1.0,0.0001);
#c1=2.0;
#c2=2.0*c1;
#c3=2.0;
#k=-1.0;
#intag(UDFF,1.0,0.0001);
for m in ak:
    k=m;
    intag(UDFF,1.0,0.0001);
    #a = np.arange(0.001,1.0,0.001);
    #plt.plot(a,UDFF(a));
plt.legend();
plt.savefig('UDF_a_vs_t.jpg',dpi=300);
plt.show();
    
for m in ak:
    k=m;
    intag(UDFz,8.0,0.0001);
    #a = np.arange(0.001,1.0,0.001);
    #plt.plot(a,UDFF(a));
plt.ylabel('Redshift z');
plt.legend();
plt.savefig('UDF_z_vs_t.jpg',dpi=300);
plt.show();

x = np. arange(0.001,15,0.1);
plt.plot(x,UDFrho(x,-0.7,-0.005,1.7005));
plt.yscale('symlog');
plt.xscale('log');
plt.ylabel(r"$\rho$");
plt.xlabel(r"a");
plt.savefig('UDF_rho.jpg',dpi=300);
plt.show();

for m in ak:
    k=m;
    x = np. arange(0.001,2,0.01);
    plt.plot(x,np.sqrt(UDFH(x)),label='k='+str(k));
plt.yscale('log');
#plt.ylim(0.8,1.1);
plt.ylabel(r"$H$");
plt.xlabel(r"a");
plt.legend();
plt.savefig('UDF_H.jpg',dpi=300);
plt.show();

x = np. arange(0.001,2,0.001);
plt.plot(x,UDFacc(x,-0.7,-0.000000000000005,1.0-0.700000000000005));
#plt.yscale('symlog');
plt.ylim(-5.0,5.0);
plt.ylabel(r"$\ddot{a}/a$");
plt.xlabel(r"a");
plt.grid(True);
plt.savefig('UDF_addot.jpg',dpi=300);
plt.show();
    
k=1.0;
x = np. arange(0.001,4,0.001);
plt.plot(x,-UDFacc(x,c3,c1,c4)/(UDFH(x)));
#plt.yscale('symlog');
plt.ylabel(r"$q$");
#plt.ylim(-0.4,0.4);
#plt.xlim(0.4,1.0);
plt.xlabel(r"a");
plt.grid(True);
plt.savefig('UDF_q.jpg',dpi=300);
plt.show();
###############################################################################
#                               Concordance                                   #
###############################################################################

print '###############################################################################\n                 Concordance model                  \n###############################################################################'; 



