#!/usr/bin/python3

# Aref Hashemi
# This python code simulates a hot object cooling down due to convection and radiation in a room
# h: convection heat transfer coefficient (W/m^2.K)
# m: mass of the object (kg)
# A: area of the body (m^2)
# Cp: specific heat capacity (J/kg.K)
# t: time (s)
# sigma: Stefan-Boltzmann constant (W/m^2.K^4)

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

g = 9.81;

# assume that the object is a metallic sphere
a = 5e-2;# radius of the sphere (m)
rho = 7800;# density of the object (kg/m^3)
k = 100;# thermal conductivity (W/m.K)
Cp = 450;

m = 4/3*math.pi*pow(a,3)*rho;
A = 4*math.pi*pow(a,2);
epsilon = 1;# set epsilon=0 for no radiation
sigma = 5.7e-8;
T_inf=298.15;# ambient temperature (K) (same as room surface temperature)
T0=T_inf+150;# initial temperature (K)
h = 20;

# tau_coduction = pow(a,2)/(k/(rho*Cp));# conduction time scale
tau_convection = m*Cp/(h*A);
# print("conduction time scale (min) = " + str(tau_coduction/60))
# print("convection time scale (min) = " + str(tau_convection/60))

mult = 3;
tf = mult*tau_convection;
dt_vec = np.array([5, 1])*60;

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, tf/60), ylim=(0, T0-T_inf))
plt.xlabel("$t\;(\mathrm{min})$")
plt.ylabel("$T-T_{\infty}\;(\mathrm{K})$")
lines = [];

for p in range(len(dt_vec)):
    dt=dt_vec[p];
    tM=math.floor(tf/dt);
    t=np.linspace(0,tf,tM);
    dt =t[1]-t[0];
    T=np.zeros(tM);
    T[0] = T0;
    for i in range(tM-1):
        T[i+1] = T[i]-dt*(h*A/(m*Cp)*(T[i]-T_inf)+sigma*epsilon*A/(m*Cp)*(pow(T[i],4)-pow(T_inf,4)));
        
    lines += plt.plot(t/60,T-T_inf, label="$" + str(round(dt/60)) + "$");


labels = [l.get_label() for l in lines]
plt.legend(lines, labels, frameon=False, title="$\Delta t\;(\mathrm{min})$")
plt.savefig('radiatingBody.pdf');
plt.show()
