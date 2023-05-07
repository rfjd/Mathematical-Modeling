#!/usr/bin/python3

# Aref Hashemi
# This python code plots the bifurcation curves for the insect outbreak problem:
# dx/dt = rx(1-x/k)-x^2/(1+x^2)

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

xM=1000;
x_min=1.01;
x_max=20  ;
x=np.linspace(x_min,x_max,xM)
k=np.zeros(xM);
r=np.zeros(xM);

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, 40), ylim=(0, 0.8))
plt.xlabel("$k$")
plt.ylabel("$r$")

for i in range(len(x)):
    r[i]=2*pow(x[i],3)/(pow(1+pow(x[i],2),2));
    k[i]=2*pow(x[i],3)/(pow(x[i],2)-1);

plt.plot(k,r);
plt.text(7,0.15, 'refuge')
plt.text(25,0.35, 'bistable')
plt.text(25,0.7, 'outbreak')
plt.savefig('insectOutbreak.pdf')
plt.show()
