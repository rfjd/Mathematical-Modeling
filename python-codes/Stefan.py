#!/usr/bin/python3

# Aref Hashemi
# This python code simulates a one-dimensional Stefan problem

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf, erfc

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

q = 336000; # latent heat of melting, J/kg
Cp = 4184; # specific heat capacity of water, J/(kg.K)
k = 0.598; # thermal conductivity of water, W/(m·K)
rho = 1000; # density of water (ice), kg/(m^3)
dT = 25; # T_W-T_0 (K)

alpha = k/(rho*Cp);
beta = Cp*dT/(q*math.sqrt(math.pi));

### calculation of lambda using Newton's method
lam_k = 1;
err = 1e10;
tol = 1e-4;
while err > tol:
    f = lam_k*erf(lam_k)*math.exp(pow(lam_k,2))-beta;
    df= erf(lam_k)*math.exp(pow(lam_k,2))+2/math.sqrt(math.pi)*lam_k+2*pow(lam_k,2)*erf(lam_k)*math.exp(pow(lam_k,2));
    lam_k += -f/df;
    err = abs(f);

lam = lam_k;
### s(t) versus time
tf = 3600 ; # seconds
tM = 10000; # number of data points to be plotted
t = np.linspace(0,tf,tM);
s = 2*lam*np.sqrt(alpha*t);

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, tf/60))
plt.xlabel("$t\;(\mathrm{min})$")
plt.ylabel("$s(t)\;(\mathrm{cm})$")

plt.plot(t/60,s*100);
plt.text(10,1, 'ice')
plt.text(10,0.25, 'water')

plt.savefig('Stefan_sVSt.pdf')
plt.show()

### theta versus x at differemt times
# vector of x values to be plotted: theta = 1-erf(x/(sqrt(4*alpha*t)))
L = 1.5*s[-1];# plotting domain (m): x \in [0,L]

xM = 1000;
x = np.linspace(0,L,xM);
# initializing the theta vector

# vector of times at which theta versus x will be plotted
t_vec = np.array([1, 10, 20, 30, 40, 50, 60])*60;# (seconds)

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, L*100), ylim=(0, 1))
plt.xlabel("$x\;(\mathrm{cm})$")
plt.ylabel("$\\theta$")

lines = [];
for p in range(len(t_vec)):
    t = t_vec[p];
    s = 2*lam*math.sqrt(alpha*t);
    index = np.where(x>s);
    index = index[0][0];
    theta = np.zeros(xM);
    theta[0:index] = 1-1/erf(lam)*erf(x[0:index]/math.sqrt(4*alpha*t));
    lines += plt.plot(x*100,theta, "-", label="$" + str(round(t/60,0)) + "$");

labels = [l.get_label() for l in lines]
plt.legend(lines, labels, frameon=False, title="$t\;(\mathrm{min})$")
plt.savefig('Stefan_thetaVSx.pdf')
plt.show()
