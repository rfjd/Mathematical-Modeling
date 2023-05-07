#!/usr/bin/python3

# Aref Hashemi
# time of evaporation for a water droplet

import math
import numpy as np

#####################################################
# params (SI Units)
rho = 1000;
MW = 18e-3;
C = 1000/MW;
Ps = 3e3;#sat pressure
T = 298.15;#K
R = 8.314;#J/mol.K
Cs = Ps/(R*T);# mol/m^3 (sat concentration)
RH = 0.5;# Relative Humidity
Cinf = Cs*RH;# mol/m^3
D = 2.5e-5;#m^2/s

a0 = 5e-6;#initial radius of the droplet (breathing, talking)
ac = 0.5*a0;#final radius of the droplet

t=(1-pow(ac/a0,2))*C*pow(a0,2)/(2*D*(Cs-Cinf));
print(str(t*1000) + " in miliseconds")
