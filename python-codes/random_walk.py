#!/usr/bin/python3

# Aref Hashemi
# This python code simulates a random walk: m steps to the right out of N total steps

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

#####################################################
# The probability of having 'm' steps to the right out of 'N' total steps  
N = 10;# total number of steps
p = 0.5;# probability of moving to the right

m = 5;

nt_max = 1000;# maximum number of trials
nt_min = 1   ;# minimum number of trials
NT = np.arange(nt_min,nt_max);# array of number of trials
P_sim = np.zeros(NT.size);# array of probabilities
for k in range(NT.size):
    count = np.zeros(NT[k], dtype=int);# will store the total number of right steps for each of the NT[k] trials
    for i in range(NT[k]):
        for j in range(N):
            s = np.random.random((1));
            if s < p:
                count[i] += 1;# one step to the right

    occ = np.count_nonzero(count == m);# Count occurrence of element 'm' in count
    P_sim[k] = occ/NT[k];# Probablity of m right steps from simulation (cumulative)

P_formula = math.factorial(N)/(math.factorial(m)*math.factorial(N-m))*pow(p, m)*pow(1-p, N-m);# Probablity of m right steps from formula

# print("Probablity of m right steps from simulation:")
# print(P_sim)
# print("Probablity of m right steps from formula:")
# print(P_formula)

Lx = 10;# in
Ly = 5 ;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(nt_min, nt_max), ylim=(0, 1))
plt.xlabel("number of trials")
plt.ylabel("$P_N(m)$")

ax.plot(NT, np.full(NT.size, P_formula))
ax.scatter(NT, P_sim, s=1, color = 'r')
ax.legend(['formula', 'simulation'], frameon=False)

plt.savefig('random_walk.pdf')
plt.show()
