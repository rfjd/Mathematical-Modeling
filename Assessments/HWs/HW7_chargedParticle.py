#!/usr/bin/python3

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

gamma = 2*math.pi;
mult = 10;
tf = mult*gamma;
dt = 0.01*gamma;
tM = math.floor(tf/dt);
t = np.linspace(0,tf,tM);
dt = t[1]-t[0];
 
epsilon = 0.05;

# regular perturbation
Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, mult))
plt.xlabel("$\\tilde{t}/2\pi$")
plt.ylabel("$\\tilde{v}$")

lines = [];

# first-order perturbation
v_perturbation=1-np.cos(t)+1/12*epsilon*(np.sin(3*t)-9*np.sin(2*t)+45*np.sin(t))-5/2*epsilon*t;

# numerical
v = np.zeros(tM);
v[0] = 0;
for i in range(tM-1):
    v[i+1] = v[i]+dt*(np.sin(t[i])-epsilon*pow(v[i],3));

lines += plt.plot(t/gamma,v,linestyle="-", color='r',label="Euler");
lines += plt.plot(t/gamma,v_perturbation,linestyle="--", color='b',label="regular perturbation");

labels = [l.get_label() for l in lines]
plt.legend(lines, labels, frameon=False)
plt.savefig('HW7_reg.pdf')
plt.show()

# two-timing
Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, mult))
plt.xlabel("$\\tilde{t}/2\pi$")
plt.ylabel("$\\tilde{v}$")

lines = [];

# two-timing first-order pertubation solution
A = np.sqrt(np.divide(3,5*np.exp(3*epsilon*t)-2));
v_twoTiming=A-np.cos(t)+epsilon*(1/12*np.sin(3*t)-3/4*A*np.sin(2*t)+3*(pow(A,2)+1/4)*np.sin(t));

lines += plt.plot(t/gamma,v,linestyle="-", color='r',label="Euler");
lines += plt.plot(t/gamma,v_twoTiming,linestyle="--", color='b',label="two-timing");

labels = [l.get_label() for l in lines]
plt.legend(lines, labels, frameon=False)
plt.savefig('HW7_twoTiming.pdf')
plt.show()
