#!/usr/bin/python3

# run by executing ./HW4_orange.py

import math
import numpy as np
import matplotlib.pyplot as plt

plt.rcParams['text.usetex'] = True
plt.rcParams.update({'font.size': 14})
 
a   = 5e-2;# radius of the orange (m)
rho = 1000;# density of the object (kg/m^3)
Cp  = 4184;# specific heat capacity (J/kg.K) 

m = 4/3*math.pi*pow(a,3)*rho;
A = 4*math.pi*pow(a,2);

epsilon = 1     ;# set epsilon=0 for no radiation
sigma   = 5.7e-8;# Stefan–Boltzmann (W/m^2.K^4)

Tinf   = 2+273.15 ;# ambient temperature (K)
T0     = 10+273.15;# ambient temperature (K)

Ts_vec = np.array([-20, -10, -5])+273.15;#effective sky temperature (K)
h = 20;# heat tranfer coefficient (W/m^2.K)

tau_convection = m*Cp/(h*A);

mult = 7;
tf = mult*tau_convection;# a multiple of tau
dt = tau_convection/100 ;# a small fraction of tau
tM = math.floor(tf/dt);

t = np.linspace(0,tf,tM);
T = np.zeros(tM);
T[0] = T0;

T_equilibrium=np.zeros(len(Ts_vec));

Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(0, tf/60))
plt.xlabel("$t\;(\mathrm{min})$")
plt.ylabel("$T\;(\mathrm{^oC})$")
lines = [];

for p in range(len(Ts_vec)):
    Ts = Ts_vec[p];

    # Newton
    T_k=Tinf;# inital guess for the equilibrium temperature
    err=1e10;
    tol=1e-7;
    while err > tol:
        f=h*(T_k-Tinf)+sigma*epsilon*(pow(T_k,4)-pow(Ts,4));
        df=h+4*sigma*epsilon*pow(T_k,3);

        T_k += -f/df;
        err = abs(f);
            
    T_equilibrium[p] = T_k;
    
    # Euler
    for i in range(tM-1):
        T[i+1] = T[i]-dt*(h*A/(m*Cp)*(T[i]-Tinf)+sigma*epsilon*A/(m*Cp)*(pow(T[i],4)-pow(Ts,4)));
        
    print("T_Newton, " + "T_Euler (Celsius)")
    print(T_equilibrium[p]-273.15,T[-1]-273.15)
    
    # lines += plt.plot(t/60,T-273.15, label="$" + str(round(Ts-273.15)) + ",\;T_{eq}=" + str(round(T_equilibrium[p]-273.15,2)) + "\;(\mathrm{^oC})$");
    lines += plt.plot(t/60,T-273.15, label="$" + str(round(Ts-273.15)) + "$");

lines += plt.plot(t/60,(Tinf-273.15)*np.ones(tM), label="$T_{\infty}$", linestyle="--", color="k");
lines += plt.plot(t/60,np.zeros(tM), label="$T_{\mathrm{freeze}}$", linestyle="-", color="k");

labels = [l.get_label() for l in lines]
plt.legend(lines, labels, frameon=False, title="$T_s\;(\mathrm{^oC})$")
plt.savefig('orange.pdf');
plt.show()


TM = 100;# number of Tinf values
Tinf_min = 0 ;
Tinf_max = 10;
Tinf_vec = np.linspace(Tinf_min,Tinf_max,TM)+273.15;
Ts_vec = np.power(h/(sigma*epsilon)*(273.15-Tinf_vec)+pow(273.15,4),0.25);
Lx = 5;# in
Ly = 5;# in
plt.figure(figsize=[Lx, Ly])
ax = plt.axes([0.15, 0.15, 0.8, 0.8], xlim=(Tinf_min, Tinf_max))
plt.xlabel("$T_{\infty}\;(\mathrm{^oC})$")
plt.ylabel("$T_{s}\;(\mathrm{^oC})$")

plt.plot(Tinf_vec-273.15,Ts_vec-273.15);
plt.savefig('Ts_VS_Tinf.pdf');
plt.show()

Tinf_Critical = 273.15+sigma*epsilon/h*pow(273.15,4);
print("critical T_inf (K) = " + str(Tinf_Critical))
