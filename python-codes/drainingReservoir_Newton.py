#!/usr/bin/python3

# Aref Hashemi
# This python code simulates a cone-shape draining reservoir using Newton's method
# h: height (head) of the liquid (m)
# r: upper radius of the liquid (m)
# a: radius of the nozzle (m)
# rho: density of liquid (kg/m^3)
# g: gravitational acceleration (m/s^2)
# dt: time step (t)
# t: time (s)
# theta: angle of the cone (rad)
# epsilon: tolerance

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

g = 9.81;
a = 1e-2;
theta = math.pi/4;
gamma = pow(a,2)*math.sqrt(2*g*math.tan(theta));
h0 = 1;
r0 = a+h0*math.tan(theta);
mult = 1.2;
tf = mult*(1/3*math.sqrt(h0/(2*g))*pow(r0/a,2));# total simulation time: a mult of ~V0/(u_nozzle*A_nozzle)
dt_vec = [60, 30, 20, 10, 5, 1];# does it matter?
epsilon = 1e-4;

C=2/15*math.sqrt(r0-a)*(8*pow(a,2)+4*a*r0+3*pow(r0,2));
r_m = 10*a;# minimum allowable radius

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, tf/60), ylim=(0, h0*100))
plt.xlabel("$t\;(\mathrm{min})$")
plt.ylabel("$h\;(\mathrm{cm})$")
lines = [];

for p in range(len(dt_vec)):
    dt=dt_vec[p];
    tM=math.floor(tf/dt);
    t=np.linspace(0,tf,tM);
    r=np.zeros(tM);
    r[0] = r0;
    h=np.zeros(tM);
    h[0] = h0;
    for i in range(tM-1):
        if r[i] < r_m:
            r = r[0:i];
            h = h[0:i];
            t = t[0:i];
            break

        r_k=r[i];
        err=1e10;
        while err > epsilon:
            f=2/15*math.sqrt(r_k-a)*(8*pow(a,2)+4*a*r_k+3*pow(r_k,2))+gamma*t[i+1]-C;
            df=2/15*(2/math.sqrt(r_k-a)*(8*pow(a,2)+4*a*r_k+3*pow(r_k,2))+math.sqrt(r_k-a)*(4*a+6*r_k));
            r_k += -f/df;
            err = abs(f);
            # print(err)
            
        r[i+1] = r_k;
        h[i+1] = (r[i+1]-a)/math.tan(theta);

    lines += plt.plot(t/60,h*100, label="$" + str(round(dt_vec[p],2)) + "$");


labels = [l.get_label() for l in lines]
plt.legend(lines, labels, frameon=False, title="$\Delta t\;(\mathrm{s})$")
plt.savefig('drainingReservoir_Newton.pdf')
plt.show()
