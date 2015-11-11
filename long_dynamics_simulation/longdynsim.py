# -*- coding: utf-8 -*-
"""
Created on Fri Oct 02 15:45:05 2015

@author: kvmu
"""

from math import pi
from math import sin
from math import cos
import numpy as np
from matplotlib import pyplot as plt

# Runge-Kutta 4th order solver, loops unraveled.
# The inputs are: the canonical conjugate varialbes phi, W, dphi_dt, dW_dt, dt.
# phi, W are floats; canonical conjugate variables
# dphi_dt, dW_dt are callable functions; can be derived from the hamiltonian
# dt is a float; time step

def RK4_step(phi, W, dphi_dt, dW_dt, dt):
    phi_1 = dphi_dt(W)*dt
    W_1 = dW_dt(phi)*dt

    phi_k = phi + phi_1*0.5
    W_k = W + W_1*0.5

    phi_2 = dphi_dt(W_k)*dt
    W_2 = dW_dt(phi_k)*dt

    phi_k = phi + phi_2*0.5
    W_k = W + W_2*0.5

    phi_3 = dphi_dt(W_k)*dt
    W_3 = dW_dt(phi_k)*dt

    phi_k = phi + phi_3
    W_k = W + W_3

    phi_4 = dphi_dt(W_k)*dt
    W_4 = dW_dt(phi_k)*dt

    phi = phi + (phi_1 + 2*(phi_2 + phi_3) + phi_4)/6
    W = W + (W_1 + 2*(W_2 + W_3) + W_4)/6
    return phi, W

def RK4_solve(phi_array, W_array, dphi_dt, dW_dt, nturns):
    frequency = 1.56e6 # Hz
    period = 1/frequency # s
    resolution_factor = 4
    dt = period/resolution_factor

    results_matrix = np.zeros((len(phi_array), nturns*resolution_factor), dtype=object)

    for i, ps_point in enumerate(zip(phi_array, W_array)):
        tempPhi, tempW = ps_point
        for step in range(nturns*resolution_factor):
                nextPhi, nextW = RK4_step(tempPhi, tempW, dphi_dt, dW_dt, dt)
                results_matrix[i, step] = (nextPhi, nextW)
                tempPhi, tempW = nextPhi, nextW

    return results_matrix

## Hamiltonian Parameters
h = 1.0 # Harmonic Number
eta = 0.86 # Dispersion
mass = 938 # MeV/c^2
R_s = 10.0 # m
q = 1.602177*10**(-19) # C
phi_s = 0.0 # phase angle
gamma = 1.010660980810234 # lorentz gamma factor
V_0 = 1000 # V
J_to_MeV = 6.241509e12 # conversion factor
c_sq = 8.988e16 # (m/s)^2

 # Constants
a = (-h*eta/(mass*R_s**2*gamma))
b = q*V_0/(2*pi)*J_to_MeV
c = sin(phi_s)
d = cos(phi_s)

# Time step
#dt = 10**(-12) # 1 ps

## Callable functions

# The Hamiltonian
def H(phi, W):
    return a/2*W**2*c_sq + b*(np.cos(phi) - d + (phi - phi_s)*c)

# Time evolution of the phase
def dphi_dt(W):
    return a*W*c_sq

# Time evolution of the angular momentum
def dW_dt(phi):
    return b*sin(phi) - b*c

# Define paramters for simulation
sqrt_num_particles = 20 # sqrt(N)

# Create the beam in phase space (it's a square right now, but doesn't matter -- physics is still the same)
phi = np.linspace(-pi/3, pi/3, sqrt_num_particles)
W = np.linspace(-1e-10, 1e-10, sqrt_num_particles)
beam = np.dstack(np.meshgrid(phi, W)).reshape(-1, 2)

# Obtain phi-W phase space
phi, W = beam[:,0], beam[:,1]

# Solve the system using 4th order Runge Kutta (4 steps = 1 turn)
results = RK4_solve(phi, W, dphi_dt, dW_dt, 1000)


#---------- PLOTTING ---------#

# Use TeX for the text on all plots
plt.rc('text', usetex=True)
plt.rc('font', family='serif')

# Generate a 100x100 grid of points for the hamiltonian contours
phi_con = np.linspace(-pi, pi, 100)
W_con = np.linspace(-1e-7, 1e-7, 100)

# Generate a mesh to be used for the contourf and contour function
phi_con, W_con = np.meshgrid(phi_con, W_con)

# Evaluate the hamiltonian for the mesh
ham = H(phi_con,W_con)

# Generate the level curves
levels = np.linspace(min(ham.flatten()),max(ham.flatten()), 100)

## Plot contours of H and visualize the phase space
# Create figure/axes object
fig, axs = plt.subplots(nrows=1, ncols=1, figsize=(10,8))

cs = axs.contourf(phi_con, W_con, ham, levels = levels[::3], cmap=plt.cm.Blues)
cs_con = axs.contour(phi_con, W_con, ham, levels=levels[::2], colors='lightgreen')

# Hardcode the x-axis to show trigonometric labels
plt.xticks([(-2 * pi), (-3 * pi/2), -pi, -pi/2, 0, pi/2, pi, (3 * pi/2), (2 * pi)],
        [r'$-2\pi$',r'$-\frac{3\pi}{2}$', r'$-\pi$',r'$-\frac{\pi}{2}$',r'$0$',
         r'$\frac{\pi}{2}$', r'$\pi$', r'$\frac{3\pi}{2}$', r'$2\pi$'], size = 15)

# Plot the colorbar (shows values of the z-axis)
colbar = fig.colorbar(cs, ax=axs, format='%.8f')
colbar.set_label('$H$ $[MeV]$ ', rotation = 270, size = 15, labelpad=30)
colbar.add_lines(cs_con)

# Synchronous condition
plt.plot(0, 0, linewidth = 4, color = 'white', marker='o' )

# Axis labels
plt.xlabel('$\phi$ $[rad]$', size= 15)
plt.ylabel('$W$ $[MeV \cdot s]$', size=15)


# This is the text for the hamiltonian, to be displayed
txtham = 'H$(q,W) = -\\frac{h \\eta}{2 m_0 R^2 \\gamma}W^2 + \\frac{qV_0}{2\\pi}(\\cos\\phi - \\cos\\phi_s + (\\phi - \\phi_s)\\sin\\phi_s)$'

# These are the parameters of the hamiltonian, to be displayed
row0 = 'Hamiltonian Paramaters:\n'
row1 = '$h=%.2f$, $\eta=%.2f$, $m_0=%.2f$ $\\frac{\\bf{MeV}}{\\bf{c^2}}$, $R=%.2f$ m\n'%(1, 0.86, 938, 10)
row2 = '$q=1.602177\\cdot10^{-19}$ C, $\phi_s=0$, $\gamma=%.2f$, $V_0=%.2f$ V '%(1.010660, 1000)
textstr = row0+row1+row2

# These are matplotlib.patch.Patch properties
props = dict(boxstyle='round', facecolor='white', alpha=0.5)

# Write the text on the plot
axs.text(0.07, 0.95, txtham, transform=axs.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

axs.text(0.1, 0.15, textstr, transform=axs.transAxes, fontsize=12,
        verticalalignment='top', bbox=props)

axs.text(0.6, -0.0725, 'Equation Source (page 21): http://trshare.triumf.ca/~craddock/TRI-87-2.pdf',
         transform=axs.transAxes, fontsize = 5, verticalalignment='top', bbox=props)

# Hiding axis ticks
plt.tick_params(axis="both", which="both", bottom="on", top="off",
                labelbottom="on", left="off", right="off", labelleft="on")

# Adding horizontal grid lines
axs.yaxis.grid(True)

# Remove axis spines
axs.spines["top"].set_visible(True)
axs.spines["right"].set_visible(False)
axs.spines["bottom"].set_visible(True)
axs.spines["left"].set_visible(False)


plt.savefig('phi-W_phasespace.pdf')
plt.show()


# Create the particle 'simulation'
# Four axes, returned as a 2-d array
f, axarr = plt.subplots(2, 2,  figsize=(12,12))

# Top left -- 0 turns // Initial beam
axarr[0, 0].scatter(phi, W)
axarr[0, 0].contour(phi_con, W_con, ham, levels=levels[::2], colors='lightgreen')
axarr[0, 0].set_ylim(-1e-8,1e-8)
axarr[0, 0].set_xlim(-pi/3-0.4,pi/3+0.4)
axarr[0, 0].set_title('0 Turns')
axarr[0, 0].set_ylabel('$W$ [$MeV \cdot s$]', size = 15)

for i in range(sqrt_num_particles**2):
#   axarr[0, 1].scatter(*zip(*results[i, 40]))
    axarr[0, 1].scatter(*results[i, 999])
axarr[0, 1].contour(phi_con, W_con, ham, levels=levels[::2], colors='lightgreen')
axarr[0, 1].set_ylim(-1e-8,1e-8)
axarr[0, 1].set_xlim(-pi/3-0.4,pi/3+0.4)
axarr[0, 1].set_title('250 Turns')

for i in range(sqrt_num_particles**2):
    axarr[1, 0].scatter(*results[i, 1999])
axarr[1, 0].contour(phi_con, W_con, ham, levels=levels[::2], colors='lightgreen')
axarr[1, 0].set_ylim(-1e-8, 1e-8)
axarr[1, 0].set_xlim(-pi/3-0.4,pi/3+0.4)
axarr[1, 0].set_title('500 Turns')
axarr[1, 0].set_ylabel('$W$ [$MeV \cdot s$]', size = 15)
axarr[1, 0].set_xlabel('$\\phi$ [$rad$]', size = 15)
#plt.sca(axarr[1,0])
#plt.xticks([-pi, -pi/2, 0, pi/2, pi],
#        [r'$-\pi$',r'$-\frac{\pi}{2}$',r'$0$',
#         r'$\frac{\pi}{2}$', r'$\pi$'], size = 15)

for i in range(sqrt_num_particles**2):
    axarr[1, 1].scatter(*results[i, 2999])
axarr[1, 1].contour(phi_con, W_con, ham, levels=levels[::2], colors='lightgreen')
axarr[1, 1].set_ylim(-1e-8,1e-8)
axarr[1, 1].set_xlim(-pi/3-0.4,pi/3+0.4)
axarr[1, 1].set_title('750 Turns')
axarr[1, 1].set_xlabel('$\\phi$ [$rad$]', size = 15)
#plt.sca(axarr[1,1])
#plt.xticks([-pi, -pi/2, 0, pi/2, pi],
#        [r'$-\pi$',r'$-\frac{\pi}{2}$',r'$0$',
#         r'$\frac{\pi}{2}$', r'$\pi$'], size = 15)

# Fine-tune figure; hide x ticks for top plots and y ticks for right plots
plt.setp([ax.get_xticklabels() for ax in axarr[0, :]], visible=False)
plt.setp([ax.get_yticklabels() for ax in axarr[:, 1]], visible=False)

plt.savefig('beam-phasespace.pdf')
plt.show()














