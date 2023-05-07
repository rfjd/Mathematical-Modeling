#!/usr/bin/python3

# Aref Hashemi
# This python code simulates the population growth problem using Euler's method
# eqn: dn/dt=rn(1-n)
# n(t): normalized population (N/K) at time t
# r: growth rate constant (1/time)
# n0: initial normalized population

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

r = 1;
mult = 6;
tf = mult*1/r;# total simulation time: a mult of ~1/r
tM = 1000;
t=np.linspace(0,tf,tM);
dt=t[1]-t[0];

numCurves=10;
n0_vec = np.linspace(0, 1.8, numCurves);
n0_vec = np.delete(n0_vec, np.where(n0_vec == 0));
n0_vec = np.delete(n0_vec, np.where(n0_vec == 1));

Lx = 6;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, tf), ylim=(0, 2))
plt.xlabel("$r\cdot t$")
plt.ylabel("$N/K$")
lines = [];
plt.plot(r*t,np.zeros(tM),'--',linewidth=2,color='k');
plt.plot(r*t,np.ones(tM),'--',linewidth=2,color='k');
for p in range(len(n0_vec)):
    n = np.zeros(tM);
    n[0] = n0_vec[p];
    for i in range(tM-1):
        n[i+1] = n[i]+dt*r*n[i]*(1-n[i]);

    plt.plot(r*t,n);
    
plt.savefig('populationGrowth.pdf')
plt.show()
