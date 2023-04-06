#!/usr/bin/python3

# Aref Hashemi
# rough calculations of the terminal velocity for a falling sphere

import math
import numpy as np

#####################################################
# params
rho_s = 1000; # density of the sphere (water droplet)
rho = 1; # density of the air
drho = rho_s-rho;
g = 9.8066;
mu = 2e-5;# viscosity of the air

D = 4e-3; # diameter of the sphere
Stokes = False;
if Stokes:
    V=drho*math.pow(D,2)*g/(18*mu);
else:
    CD = 0.4;
    V=math.sqrt(4*drho*D*g/(3*CD*rho));

print(V)
Re = rho*V*D/mu
print(Re)
