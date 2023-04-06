#!/usr/bin/python3

# Aref Hashemi
# This python code simulates the movement of a projcetile launched upward
# dv/dt=-g-delta*v^3
# v: velocity (m/s)
# v0: initial velocity (m/s)
# z: height (m)
# t: time (s)

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

g = 9.81;
v0 = 10;
mult = 2;
tf = mult*v0/g;
tM = 1000;
dt = tf/tM;
delta_vec = [0, 0.01, 0.1];# epsilon/m

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, tf), ylim=(0, 0.51*pow(v0,2)/g))
plt.xlabel("$t\;(\mathrm{s})$")
plt.ylabel("$z\;(\mathrm{m})$")

t=np.linspace(0,tf,tM);
z_analytical=-0.5*g*np.power(t,2)+v0*t;

# print("-----------------")
# print(t)
# print(v0*t)
# print("-----------------")

lines = [];
for p in range(len(delta_vec)):
    delta=delta_vec[p];
    v=np.zeros(tM);
    v[0] = v0;
    z=np.zeros(tM);
    z[0] = 0;
    for i in range(tM-1):
        v[i+1] = v[i]-dt*(g+delta*pow(v[i],3));
        z[i+1] = z[i]+dt*v[i];

    lines += plt.plot(t,z, "-", label="$" + str(round(delta,2)) + "$");

plt.plot(t,z_analytical, "--");

labels = [l.get_label() for l in lines]
plt.legend(lines, labels, frameon=False, title="$\delta$")
plt.savefig('projectile.pdf')
plt.show()
