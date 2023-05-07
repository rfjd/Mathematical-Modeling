#!/usr/bin/python3

import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

#### ab>0
a = 1;
b = 1;
R0_vec=np.array([ 11, 13,-11,-13,-10,-11, 10, 11]);
J0_vec=np.array([-10,-11, 10, 11, 11, 13,-11,-13]);

tf = 10;
tM = 1000;
t = np.linspace(0,tf,tM);
R = np.zeros(tM);
J = np.zeros(tM);

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=[-10,10], ylim=[-10,10])
plt.xlabel("$R$")
plt.ylabel("$J$")
lines = [];

R_plot=np.linspace(-10,10,100);
lines+=plt.plot(R_plot,math.sqrt(b/a)*R_plot, 'k', label="$\lambda_1\!=\!+\sqrt{ab}$")
lines+=plt.plot(R_plot,-math.sqrt(b/a)*R_plot, 'k', label="$\lambda_2\!=\!-\sqrt{ab}$")
R_target = 2;
dR = 1;
plt.arrow(R_target, math.sqrt(b/a)*R_target, dR, math.sqrt(b/a)*dR, shape='full', color='k', lw=1, length_includes_head=True, head_width=0.5)
plt.arrow(-R_target, -math.sqrt(b/a)*R_target, -dR, -math.sqrt(b/a)*dR, shape='full', color='k', lw=1, length_includes_head=True, head_width=0.5)
R_target+=dR;
plt.arrow(R_target, -math.sqrt(b/a)*R_target, -dR, math.sqrt(b/a)*dR, shape='full', color='k', lw=1, length_includes_head=True, head_width=0.5)
plt.arrow(-R_target, math.sqrt(b/a)*R_target, dR, -math.sqrt(b/a)*dR, shape='full', color='k', lw=1, length_includes_head=True, head_width=0.5)

A = np.array([[1, 1], [math.sqrt(b/a),-math.sqrt(b/a)]])

for p in range(len(R0_vec)):
    R0 = R0_vec[p];
    J0 = J0_vec[p];

    B = np.array([R0, J0])
    c = np.linalg.solve(A, B)
    
    R=c[0]*np.exp(math.sqrt(a*b)*t)+c[1]*np.exp(-math.sqrt(a*b)*t);
    J=c[0]*np.exp(math.sqrt(a*b)*t)*math.sqrt(b/a)-c[1]*np.exp(-math.sqrt(a*b)*t)*math.sqrt(b/a);
    
    R = np.real(R);
    J = np.real(J);

    plt.plot(R,J,'b');

    i=100
    s=1;
    plt.arrow(R[i], J[i], R[i+s]-R[i], J[i+s]-J[i], shape='full', color='b', lw=1, length_includes_head=True, head_width=0.5)


plt.scatter(0, 0, s=80, facecolor='none', edgecolors='k')
# labels = [l.get_label() for l in lines]
# plt.legend(lines, labels, frameon=False, title="eigen directions")
plt.xticks([0],labels=['$0$'])
plt.yticks([0],labels=['$0$'])
plt.savefig('HW6_ii_plot_saddle.pdf');
plt.show()


#### ab<0
a = 1;
b = -1;
R0_vec=np.array([2, 4, 6, 8]);
J0_vec=np.array([0, 0, 0, 0]);

tf = 10;
tM = 1000;
t = np.linspace(0,tf,tM);
R = np.zeros(tM);
J = np.zeros(tM);

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=[-9,9], ylim=[-9,9])
plt.xlabel("$R$")
plt.ylabel("$J$")
lines = [];

A = np.array([[1, 1], [cmath.sqrt(b/a),-cmath.sqrt(b/a)]])

for p in range(len(R0_vec)):
    R0 = R0_vec[p];
    J0 = J0_vec[p];

    B = np.array([R0, J0])
    c = np.linalg.solve(A, B)
    
    R=c[0]*np.exp(cmath.sqrt(a*b)*t)+c[1]*np.exp(-cmath.sqrt(a*b)*t);
    J=c[0]*np.exp(cmath.sqrt(a*b)*t)*cmath.sqrt(b/a)-c[1]*np.exp(-cmath.sqrt(a*b)*t)*cmath.sqrt(b/a);
    
    R = np.real(R);
    J = np.real(J);

    plt.plot(R,J,'b');

    i=0
    s=1;
    plt.arrow(R[i], J[i], R[i+s]-R[i], J[i+s]-J[i], shape='full', color='b', lw=1, length_includes_head=True, head_width=0.35)

plt.scatter(0, 0, s=80, facecolor='k', edgecolors='k')
# labels = [l.get_label() for l in lines]
# plt.legend(lines, labels, frameon=False, title="eigen directions")
plt.xticks([0],labels=['$0$'])
plt.yticks([0],labels=['$0$'])
plt.savefig('HW6_ii_plot_center.pdf');
plt.show()

