#!/usr/bin/python3

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})

R0_vec=np.linspace(0.0001,0.001,10);
R0_vec=R0_vec[::-1];
J0_vec=np.zeros(len(R0_vec));

tf = 50;
tM = 1000;
t = np.linspace(0,tf,tM);
R = np.zeros(tM);
J = np.zeros(tM);

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(-50, 50), ylim=[-50, 100])
plt.xlabel("$R$")
plt.ylabel("$J$")
lines = [];

A = np.array([[1, 1], [0.5+(math.sqrt(3)/2)*1j, 0.5-+(math.sqrt(3)/2)*1j]])

for p in range(len(R0_vec)):
    R0 = R0_vec[p];
    J0 = J0_vec[p];

    B = np.array([R0, J0])
    c = np.linalg.solve(A, B)

    R=c[0]*np.exp((0.5+(math.sqrt(3)/2)*1j)*t)+c[1]*np.exp((0.5-(math.sqrt(3)/2)*1j)*t);
    J=c[0]*np.exp((0.5+(math.sqrt(3)/2)*1j)*t)*(0.5+(math.sqrt(3)/2)*1j)+c[1]*np.exp((0.5-(math.sqrt(3)/2)*1j)*t)*(0.5-(math.sqrt(3)/2)*1j);
    R = np.real(R);
    J = np.real(J);

    lines += plt.plot(R,J,'b');
    i=445;
    s=1;
    plt.arrow(R[i], J[i], R[i+s]-R[i], J[i+s]-J[i], shape='full', color='b', lw=1, length_includes_head=True, head_width=2.5)

# R_target=20;
# i=np.argmax(R>R_target);
# s=1;
# plt.arrow(R[i], J[i], R[i+s]-R[i], J[i+s]-J[i], shape='full', color='b', lw=1, length_includes_head=True, head_width=2.5)
# plt.scatter(0, 0, s=80, facecolor='none', edgecolors='k')
plt.xticks([0],labels=['$0$'])
plt.yticks([0],labels=['$0$'])
plt.savefig('HW6_i_plot.pdf');
plt.show()
