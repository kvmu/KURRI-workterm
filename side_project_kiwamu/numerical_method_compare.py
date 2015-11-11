# -*- coding: utf-8 -*-
"""
Created on Wed Sep 02 14:02:35 2015

@author: Kevin
"""

import numpy as np
from matplotlib import pyplot as plt

E_mat = np.genfromtxt('lifetime135Cs2 .txt')
RK4_mat = np.genfromtxt('rk4solution.txt')

time_axis = RK4_mat[:,0] # E_mat and RK4_mat have the same time_axis

########################### Check 1: conservation of mass

e_cons = np.sum(E_mat[:,1:], axis=1)
rk4_cons = np.sum(RK4_mat[:,1:], axis=1)

plt.figure(figsize=(8,6))
ax = plt.subplot(111)

plt.plot(time_axis, np.abs(1-e_cons), '-',
         color='lightblue', alpha = 0.6,
         linewidth = 2, label = 'Euler Mass Cons. Error')


plt.plot(time_axis, np.abs(1-rk4_cons), '-',
         color='green', alpha = 0.6,
         linewidth = 2, label = 'Runge-Kutta Mass Cons. Error')

# labels
plt.xlabel('Time [days]', fontsize = 14)
plt.ylabel('|$1 - $total amount of Moles|', fontsize = 14)
plt.legend(loc='best', fancybox=True, fontsize=10)
plt.title('Euler (FORTRAN) v.s. Runge-Kutta (Python 2.7)')
 # hiding axis ticks
plt.tick_params(axis="both", which="both", bottom="off", top="off",
        labelbottom="on", left="off", right="off", labelleft="on")

# adding horizontal grid lines
ax.yaxis.grid(True)
ax.set_yscale('log')

# remove axis spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)

plt.savefig("cons_mass_error_compare.pdf")
plt.show()

########################### Check 2: relative error

plt.figure(figsize=(8,6))
ax = plt.subplot(111)

plt.plot(time_axis, np.abs(E_mat[:,1] - RK4_mat[:,1]), '-',
         alpha = 0.6, linewidth = 2, label = 'Cs-135')

plt.plot(time_axis, np.abs(E_mat[:,2] - RK4_mat[:,2]), '-',
         alpha = 0.6, linewidth = 2, label = 'Cs-133')

plt.plot(time_axis, np.abs(E_mat[:,3] - RK4_mat[:,3]), '-',
         alpha = 0.6, linewidth = 2, label = 'Xe-135')

plt.plot(time_axis, np.abs(E_mat[:,4] - RK4_mat[:,4]), '-',
         alpha = 0.6, linewidth = 2, label = 'Xe-134')

plt.plot(time_axis, np.abs(E_mat[:,5] - RK4_mat[:,5]), '-',
         alpha = 0.6, linewidth = 2, label = 'Xe-133')

plt.plot(time_axis, np.abs(E_mat[:,6] - RK4_mat[:,6]), '-',
         alpha = 0.6, linewidth = 2, label = 'Xe-132')

plt.plot(time_axis, np.abs(E_mat[:,7] - RK4_mat[:,7]), '-',
         alpha = 0.6, linewidth = 2, label = 'Ba-135')

# labels
plt.xlabel('Time [days]', fontsize = 14)
plt.ylabel('|Euler - Runge-Kutta|', fontsize = 14)
plt.legend(loc='best', fancybox=True, fontsize=10)
plt.title('Euler (FORTRAN) v.s. Runge-Kutta (Python 2.7)', fontsize = 16)
 # hiding axis ticks
plt.tick_params(axis="both", which="both", bottom="off", top="off",
        labelbottom="on", left="off", right="off", labelleft="on")

# adding horizontal grid lines
ax.yaxis.grid(True)
ax.set_yscale('log')

# remove axis spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)

plt.savefig("relative_error_compare.pdf")
plt.show()


########################### Check 3: neutron production relative error

# Below will be mols/sec for that day
neutron_dataRK4 = np.genfromtxt('neutron_production_rate.txt')[:,1]*60
total_neutronRK4 = np.cumsum(neutron_dataRK4*1440)

total_neutronE = np.genfromtxt('lifetime135Cs-neu2.txt')

plt.figure(figsize=(8,6))
ax = plt.subplot(111)

plt.plot(time_axis, np.abs(total_neutronE[:,1] - total_neutronRK4), '-',
         alpha = 0.6, linewidth = 2, label = 'Relative Error')

# labels
plt.xlabel('Time [days]', fontsize = 14)
plt.ylabel('|Euler - Runge-Kutta|', fontsize = 14)
plt.legend(loc='best', fancybox=True, fontsize=10)
plt.title('Euler (FORTRAN) v.s. Runge-Kutta (Python 2.7)', fontsize = 16)
 # hiding axis ticks
plt.tick_params(axis="both", which="both", bottom="off", top="off",
        labelbottom="on", left="off", right="off", labelleft="on")

# adding horizontal grid lines
ax.yaxis.grid(True)
ax.set_yscale('log')

# remove axis spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)

plt.savefig("relative_error_compare_neutron.pdf")
plt.show()


########################### Result 1: Cs-135 mol vs. Total Neutron



plt.figure(figsize=(8,6))
ax = plt.subplot(111)

plt.plot( total_neutronE[:,1], E_mat[:,1] , '-',
         alpha = 0.6, linewidth = 2, label = 'y = -0.6427605x + 0.99984943')

# labels
plt.xlabel('Total Neutrons Produced [mol]', fontsize = 14)
plt.ylabel('Cs-135 [mol]', fontsize = 14)
plt.legend(loc='best', fancybox=True, fontsize=10)
plt.title('Euler (FORTRAN)', fontsize = 16)
 # hiding axis ticks
plt.tick_params(axis="both", which="both", bottom="off", top="off",
        labelbottom="on", left="off", right="off", labelleft="on")

# adding horizontal grid lines
ax.yaxis.grid(True)

# remove axis spines
ax.spines["top"].set_visible(False)
ax.spines["right"].set_visible(False)
ax.spines["bottom"].set_visible(False)
ax.spines["left"].set_visible(False)

plt.savefig("cs_vs_neu.pdf")
plt.show()
