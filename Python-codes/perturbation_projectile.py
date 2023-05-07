#!/usr/bin/python3

# Aref Hashemi
# This python code simulates the movement of a projcetile launched upward: dv/dt=-1-eps*v^3
# v: dimensionless velocity (-) scaled by the initial velocity v0
# t: dimensionless time (-) scaled by v0/g

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

tf = 5;
tM = 1000;
dt = tf/tM;
eps_vec = [0.001, 0.01, 0.1];

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, tf))
plt.xlabel("$\\tilde{t}$")
plt.ylabel("$\\tilde{v}$")

t = np.linspace(0,tf,tM);

colors=['r','b','g','k'];
lines = [];
for p in range(len(eps_vec)):
    epsilon = eps_vec[p];

    # First-order perturbation
    v_perturbation=1-t+0.25*epsilon*(np.power(t-1,4)-1);
    # Numerical
    v = np.zeros(tM);
    v[0] = 1;
    for i in range(tM-1):
        v[i+1] = v[i]-dt*(1+epsilon*pow(v[i],3));

    lines += plt.plot(t,v,linestyle="-", color=colors[p],label="$" + str(round(epsilon,3)) + "$");
    plt.plot(t,v_perturbation,linestyle="--", color=colors[p]);


labels = [l.get_label() for l in lines]
plt.legend(lines, labels, frameon=False, title="$\epsilon$")
plt.savefig('perturbation_projectile.pdf')
plt.show()
