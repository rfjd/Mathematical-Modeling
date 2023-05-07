#!/usr/bin/python3

# run by executing ./HW5_popGrowth.py

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 22})


# (i)
a = 1;
mult = 6;
tf = mult*1/a;# total simulation time: a mult of ~1/r
tM = 1000;
t=np.linspace(0,tf,tM);
dt=t[1]-t[0];
b = 1/4;

numCurves=13;
N0_vec = np.linspace(0, 6, numCurves);
N0_vec = np.delete(N0_vec, np.where(N0_vec == 0));
N0_vec = np.delete(N0_vec, np.where(N0_vec == 1/b));

Lx = 6;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, tf), ylim=(0, 6))
plt.xlabel("$t$")
plt.ylabel("$N$")
lines = [];
plt.plot(a*t,np.zeros(tM),'--',linewidth=2,color='k');
plt.plot(a*t,1/b*np.ones(tM),'--',linewidth=2,color='k');
for p in range(len(N0_vec)):
    N = np.zeros(tM);
    N[0] = N0_vec[p];
    for i in range(tM-1):
        N[i+1] = N[i]-dt*a*N[i]*math.log(b*N[i]);

    plt.plot(a*t,N);

plt.xticks([0])
plt.yticks([0,1/b],labels=['$0$','$1/b$'])
plt.savefig('popGrowth_i.pdf')
plt.show()

# (ii)
r=4;
a=1;
b=math.sqrt(r/a);

mult = 15;
tf = mult*1/r;# total simulation time: a mult of ~1/r
tM = 1000;
t=np.linspace(0,tf,tM);
dt=t[1]-t[0];

numCurves=13;
N0_vec = np.linspace(0, 6, numCurves);
N0_vec = np.delete(N0_vec, np.where(N0_vec == 0));
N0_vec = np.delete(N0_vec, np.where(N0_vec == 2*b));

Lx = 6;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, tf), ylim=(0, 6))
plt.xlabel("$t$")
plt.ylabel("$N$")
lines = [];
plt.plot(a*t,np.zeros(tM),'--',linewidth=2,color='k');
plt.plot(a*t,2*b*np.ones(tM),'--',linewidth=2,color='k');
for p in range(len(N0_vec)):
    N = np.zeros(tM);
    N[0] = N0_vec[p];
    for i in range(tM-1):
        N[i+1] = N[i]+dt*N[i]*(r-a*pow((N[i]-b),2));

    plt.plot(r*t,N);

plt.xticks([0])
plt.yticks([0,2*b],labels=['$0$','$2b$'])
plt.savefig('popGrowth_ii.pdf')
plt.show()
