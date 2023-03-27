#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computational physics project 1
The code to calculate the pair correlation, average pressure and to plot the 
energy and pair correlation
Version 24-03-2022
"""

import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation
filename = '22_03_24_solid'

rho = 0.5
sim_time = 100

h = 0.01
N_particles = 108
dim = 3
L_sys = (N_particles/rho)**(1/3)
epsilon = 119.8 #characteristic energy for Argon in K

bins = np.loadtxt('bins'+filename+'.txt', delimiter = ',')
hist = np.loadtxt('histogram'+filename+'.txt', delimiter = ',')
hist_list = np.array([hist])
#in case you want to average over multiple simulations: load in the additional
#files and add the data to hist_list

#the pressures of 10 simulations
press_list = np.array([1.0112214555390266, 1.0110052814886874, 1.000604712212472, 1.108481078321158, 0.9995061197944768, 1.0232418707497846, 0.9098521074779538, 1.0007673170606597, 0.9136829843246153, 0.9517080439770439])

kin = np.loadtxt('kineticenergy'+filename+'.txt', delimiter =',')
LJ_pot = np.loadtxt('LJ_potential'+filename+'.txt', delimiter =',')
elec_pot = np.loadtxt('elec_potential'+filename+'.txt', delimiter = ',')
time = np.arange(0, len(kin))*h #array of time steps

traj = np.loadtxt('traj'+filename+'.txt') #Loads the data
traj = traj.reshape((N_particles, sim_time, dim)) #Reshapes it into the usual form

def paircorrelation(hist, bins):
    """"To calculate the pair correlation.
    Args:
        hist (array): the histograms of multiple simulations.
        bins (array): the wanted bins for the histogram.
    Returns:
        r: the center of the bins
        n_ave: the average of the histograms.
        g: the pair correlation function."""
    
    r = np.zeros(len(bins)-1)
    dr = np.abs(bins[0]-bins[1])
    n_ave = np.average(hist, axis = 0)
    V = L_sys**3
    g = np.zeros(len(bins)-1)
    for i in range(0, len(n_ave)):
        r[i] = (bins[i+1] + bins[i])/2
        g[i] = V*n_ave[i]/(N_particles*(N_particles-1)*4*np.pi*r[i]**2*dr)
    
    return r, n_ave, g

bins_plot, n_ave, g = paircorrelation(hist_list, bins)

#the energy plot
plt.figure()
plt.plot(time, kin, label = 'kinetic energy')
plt.plot(time, LJ_pot, label = 'lennard-jones potential energy')
plt.plot(time, elec_pot, label = 'electrical potential energy')
plt.plot(time, kin+LJ_pot+elec_pot, label = 'total energy')
plt.xlabel('time step (a.u.)')
plt.ylabel('energy ($\epsilon$)')
plt.legend()
plt.title('system energy')
plt.savefig('energy'+filename)

#particle number histogram
plt.figure()
plt.bar(bins_plot, n_ave, width = 0.1, align = 'center')
plt.xlabel('Distance between particles')
plt.ylabel('Number of particles')
plt.title('Histogram of the particle number')
plt.savefig('hist'+filename)

plt.figure()
plt.scatter(bins_plot, g)
plt.plot(bins_plot, g, '--')
plt.xlabel('Distance between particles ($\sigma$)')
plt.ylabel('g(r)')
plt.title('pair correlation')
plt.savefig('correlation'+filename)

fig = plt.figure()
ax = fig.add_subplot(111,projection='3d')
sc = ax.scatter([],[],[], c='darkblue', alpha=0.5)

def update(i):
    sc._offsets3d = (np.ndarray.tolist(traj[:,i,0]), np.ndarray.tolist(traj[:,i,1]), np.ndarray.tolist(traj[:,i,2]))

ax.set_xlabel('X')
ax.set_ylabel('Y')
ax.set_zlabel('Z')
ax.set_xlim(-0,np.amax(traj))
ax.set_ylim(-0,np.amax(traj))
ax.set_zlim(-0,np.amax(traj))

ani = matplotlib.animation.FuncAnimation(fig, update, frames=sim_time, interval=1)

plt.tight_layout()
plt.show()