import numpy as np
import matplotlib.pyplot as plt

N = 10000000
a = np.linspace(0.001,3,N)
y = np.ones(shape = (N,1))
A = 0.25
B = 0.2
C = 0.15
rho_d = A*(1/a**3)
rho_r = B*(1/a**4)
rho_lambda = C*y

fig = plt.figure()
f, ax = plt.subplots(1)
plt.xlabel('scale factor a(t)')
plt.ylabel('Density rho')
plt.title('Densities vs the scale factor with random constants')
plt.plot(a, rho_d,label="Matter C=2.5")
plt.plot(a, rho_r,label="Radiation C=2")
plt.plot(a, rho_lambda,label="Dark Energy C=1.5")
ax.set_ylim(ymin=0, ymax=1.5)
ax.set_xlim(xmin=0.5, xmax=1.5)
plt.grid(True)
plt.legend(bbox_to_anchor=(0.95, 0.95), loc=1, borderaxespad=0.)
plt.savefig("Densities vs scale factor.pdf")
plt.show(fig)