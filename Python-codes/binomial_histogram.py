#!/usr/bin/python3

# Aref Hashemi
# This python code uses the concept of random walk to find a probability distribution function

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True

#####################################################
# The probability of having 'm' steps to the right out of 'N' total steps  
N = 100;# total number of steps
p = 0.7;# probability of moving to the right


m_vec = np.arange(0,N+1);# vector of possible number of right steps

nt = 3000;# total number of trials
P_sim = np.zeros(m_vec.size);# array of probabilities from simulation
P_formula = np.zeros(m_vec.size);# array of probabilities from formula
for k in range(m_vec.size):
    m = m_vec[k]
    P_formula[k] = math.factorial(N)/(math.factorial(m)*math.factorial(N-m))*pow(p, m)*pow(1-p, N-m);# Probablity of m right steps from formula
    count = np.zeros(nt, dtype=int);# will store the total number of right steps for each of the nt trials
    for i in range(nt):
        for j in range(N):
            s = np.random.random((1));
            if s < p:
                count[i] += 1;# one step to the right

    occ = np.count_nonzero(count == m);# Count occurrence of element 'm' in count
    P_sim[k] = occ/nt;# Probablity of m right steps from simulation

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, N), ylim=(0, 0.3))
plt.xlabel("$m$")
plt.ylabel("$PDF$")

ax.plot(m_vec, P_formula)
ax.scatter(m_vec, P_sim, s=2, color = 'r')
ax.legend(['formula', 'simulation'], frameon=False)

# area under the curve
sum_sim = 0;
sum_formula = 0;
for i in range(m_vec.size):
    sum_sim += P_sim[i]
    sum_formula += P_formula[i]

print([sum_sim, sum_formula])

plt.savefig('binomial_histogram.pdf')
plt.show()
