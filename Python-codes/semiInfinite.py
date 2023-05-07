#!/usr/bin/python3

# Aref Hashemi
# This python code plots the temperature versus time for a semi-infinite domain
# theta = 1-erf(eta); theta = T-T0/(TW-T0); eta = x/sqrt(4*alpha*t)
# T: temperature
# T0: initial temperature
# TW: wall temperature
# x: distance from the wall
# alpha: thermal diffusivity
# t: time

import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import erf, erfc

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

### Dimensionless plot
# vector of eta values to be plotted: theta = 1-erf(eta)
eta_max = 2;
etaM = 100;
eta = np.linspace(0,eta_max,etaM)

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, eta_max), ylim=(0, 1))
plt.xlabel("$\\eta$")
plt.ylabel("$\\theta$")

theta = 1-erf(eta);
plt.plot(eta,theta, "-");

plt.savefig('semiInfinite_dimensionless.pdf')
plt.show()

### Dimensional plot
# vector of x values to be plotted: theta = 1-erf(x/(sqrt(4*alpha*t)))
L = 1;# plotting domain (m): x \in [0,L]
xM = 100;
x = np.linspace(0,L,xM);
# initializing the theta vector
theta = np.zeros(xM);

alpha = 2.3e-5;# thermal diffusivity of iron (m^2/s)
# vector of times at which theta versus x will be plotted
tau = pow(L,2)/alpha;
# print("time scale (min) = " + str(tau/60))
t_vec = np.linspace(1,20,10)*60;# (seconds)

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, L*100), ylim=(0, 1))
plt.xlabel("$x\;(\mathrm{cm})$")
plt.ylabel("$\\theta$")

lines = [];
for p in range(len(t_vec)):
    t = t_vec[p];
    theta = 1-erf(x/math.sqrt(4*alpha*t));
    
    lines += plt.plot(x*100,theta, "-", label="$" + str(round(t/60,0)) + "$");

labels = [l.get_label() for l in lines]
plt.legend(lines, labels, frameon=False, title="$t\;(\mathrm{min})$")
plt.savefig('semiInfinite_dimensional.pdf')
plt.show()
