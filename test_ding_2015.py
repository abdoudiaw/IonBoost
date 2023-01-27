import numpy  as np


def electron_Debye_length(electron_density,electron_temperature):
    
    return 7.43e2*np.sqrt(electron_temperature/electron_density) # cm



# Simulations parameters:
electron_temperature = Ti=31.0 # eV
electron_density =  1.5e13 # 1/cc  --> this density seems quite hight (1.5e19 1/cc)?


L_sim = 3.0 # cm
L_plasma = 1.0 # cm
L_sheath_half = 1.0  # system is symmetric overall L_sheath=2.cm

#Simulation size in reduced units:
Debye_Length= electron_Debye_length(electron_density,electron_temperature)
print("Debye_Length [cm]",Debye_Length)


L_reduced  =  L_sim/Debye_Length

L_plasma  =  L_plasma/Debye_Length
print("Computational domain in reduced units",L_reduced)

print("L_plasma",L_plasma)

q=  1.6022e-19
kB=1.38e-23

E0 = kB*electron_temperature*11600/q
print("E0",E0)
