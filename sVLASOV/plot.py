import numpy as np

import matplotlib.pyplot as plt


data = np.loadtxt('profil', skiprows=1, unpack=True) #x0      xt      v       phi     E       ni      ne      rho     Te      vemoy

# Plot density profile
plt.plot(data[0], data[5], label='ion density')
plt.plot(data[0], data[6], label='electron density')
plt.xlim(-20,20)
plt.xlabel('x [m]')
plt.ylabel('density [m^-3]')
plt.legend()
plt.show()

# read distribution function

fe = np.loadtxt('fe', skiprows=1, unpack=True) # time    ve(j)   fe0(j)  E(j)    Te(j)
# plot distribution function as function of E
plt.plot(fe[3], fe[2])
plt.xlabel('E [eV]')
plt.ylabel('fe0')
plt.show()
