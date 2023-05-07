#!/usr/bin/python3

# Aref Hashemi
# This python code simulates weakly damped oscillator: d^2x/dt^2+2*eps*dx/dt+x=0

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

tf = 50;
tM = 1000;
epsilon = 0.1;

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, tf))
plt.xlabel("$t$")
plt.ylabel("$x$")

t = np.linspace(0,tf,tM);

twoTiming=True;

# First-order perturbation
x_perturbation=np.sin(t)-epsilon*np.multiply(t,np.sin(t));

# Two-timing
x_twoTiming=np.multiply(np.exp(-epsilon*t),np.sin(t));

# Analytical
omega=pow(1-pow(epsilon,2),1/2);
x_analytical=np.multiply(np.exp(-epsilon*t),np.sin(omega*t))/omega;


lines = [];
lines += plt.plot(t,x_analytical,"-r",label="analytical");
lines += plt.plot(t,x_perturbation,"-b",label="perturbation");
if twoTiming:
    lines += plt.plot(t,x_twoTiming,"--g",label="two-timing");


labels = [l.get_label() for l in lines]
plt.legend(lines, labels, frameon=False)
plt.savefig('perturbation_damped.pdf')
plt.show()
