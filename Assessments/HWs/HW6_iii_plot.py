#!/usr/bin/python3

import math
import cmath
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

#### ab>0
a = -1;
b = 1;
R0_vec=np.array([ 10, 10, 10, 10, 10,  8,  4,  0, -4, -8,-10,-10,-10,-10,-10,  8,  4,  0, -4, -8]);
J0_vec=np.array([  8,  4,  0, -4, -8, 10, 10, 10, 10, 10,  8,  4,  0, -4, -8,-10,-10,-10,-10,-10]);

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
lines+=plt.plot(R_plot,R_plot, 'k', label="$\lambda_1\!=\!+\sqrt{ab}$")
lines+=plt.plot(R_plot,-R_plot, 'k', label="$\lambda_2\!=\!-\sqrt{ab}$")
R_target = 3;
dR = 1;
R_target+=dR;
plt.arrow(R_target, -R_target, -dR, dR, shape='full', color='k', lw=1, length_includes_head=True, head_width=0.5)
plt.arrow(-R_target, R_target, dR, -dR, shape='full', color='k', lw=1, length_includes_head=True, head_width=0.5)

A = np.array([[1, 1], [1,-1]])

for p in range(len(R0_vec)):
    R0 = R0_vec[p];
    J0 = J0_vec[p];

    B = np.array([R0, J0])
    c = np.linalg.solve(A, B)
    
    R=c[0]+c[1]*np.exp((a-b)*t);
    J=c[0]-c[1]*np.exp((a-b)*t);
    
    R = np.real(R);
    J = np.real(J);

    plt.plot(R,J,'b');
    point=0.5*(R0+J0);
    plt.scatter(point, point, s=80, facecolor='k', edgecolors='k')

    i=100
    s=1;
    plt.arrow(R[i], J[i], R[i+s]-R[i], J[i+s]-J[i], shape='full', color='b', lw=1, length_includes_head=True, head_width=0.5)

plt.scatter(0, 0, s=80, facecolor='k', edgecolors='k')
# labels = [l.get_label() for l in lines]
# plt.legend(lines, labels, frameon=False, title="eigen directions")
plt.xticks([0],labels=['$0$'])
plt.yticks([0],labels=['$0$'])
plt.savefig('HW6_iii_plot.pdf');
plt.show()
