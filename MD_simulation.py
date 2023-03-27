#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Computational physics project 1
The code to run the simulation and save the data
Version 24-03-2022
"""

import numpy as np

#the simulation constants:
Temp = 1.2 #temperature (in K)
rho = 0.5 #density
filename = '22_03_24_solid'+'.txt'

sim_time = 100 #number of time steps in the simulation
E_field = 0.0 #electrical field 
E_charge = 0 #electric charge
#Try E_field = 5.0, E_charge = 1 in the solid phase...

N_particles = 108 
Dim = 3
L_sys = (N_particles/rho)**(1/3) #size of the system
test_time = 10 #number of time steps after which to evaluate the kinetic energy and scale factor
threshold = 0.01 #allowed deviation from the scale factor
h = 0.01 #time step
k_c = 1 / (4 * np.pi * 1) #constant in front of Coulomb Force
sigma = np.sqrt(Temp) #width for the gaussian distribution
bins = np.arange(0, L_sys/2, 0.1) #bins for the histogram of the distances between the particles

def moveparticles(pos_i, vel_i, E_charge):
    """"
    To calculate the new position of the particles after one 
    time step using the velocity-verlet (vv) algorithm.
    
    Args:
        pos_i (array): Positions of each particle in 
        compoments at this timestep. 
        vel_i (array): Velocities of each particle in components
        at this timestep.
        E_charge (integer): The particles electric charge.
        
    Returns:
        pos%L_sys (array): Positions of each particle in 
        compoments at the new timestep. Division modulo (%) 
        L_sys implements the periodic boundary condistions. 
        vel (array): Velocities of each particle in components
        at the new timestep. 
        r (array): The distances between the particles.
        pres (integer): Parameter used to calculate pressure.
    """
    F_i, r, pres = particleforce(pos_i, E_charge)
    pos_new = pos_i + vel_i*h + h**2*np.sum(F_i, axis = 1)/2
    F_new, r, pres = particleforce(pos_new, E_charge)
    vel_new = velocities(vel_i, F_i, F_new)
    return pos_new%L_sys, vel_new, r, pres

def velocities(vel_i, F_i, F_new): #version 2
    """To calculate the velocities with the vv algorithm.
    Args:
        pos (array): The position at the current timestep. #i think we do not actually need this parameter in this function
        v0 (array): The velocity at the current timestep.
        F (array): The force of the previous timestep.
        F_new (array): The force of the current timestep.
    
    Returns:
        vnew (array): The new velocity.
    """
    vel_new = vel_i
    for i in range(N_particles):
      vel_new[i] = vel_i[i] + (np.sum(F_new[i], axis = 0) + np.sum(F_i[i], axis = 0))*h/2
    return vel_new

def particleforce(pos, E_charge):
    """To calculate the Lennard-Jones force between the particles.
    Args: 
        pos (array): The position at the current timestep.
        E_field (integer): The magnitude of the electric field.
    Returns:
        F (array): The force between each pair of particles.
        r_array (array): The distance between each pair of particles.
        pres (scalar): The pressure of the system (/2 corrects for double counting)
        """
    F = np.zeros((N_particles, N_particles-1, Dim))
    r_array = np.zeros((N_particles, N_particles-1)) #for the histogram
    pres = 0 #for calculating the pressure
    for i in range(N_particles):
      z = 0 # a counter
      for j in [x for x in range(N_particles) if x != i]: #list comprehension, excludes the one particle we are calculating the force "from"
        r = (pos[i] - pos[j] + L_sys/2)%L_sys - L_sys/2
        F[i][z] = LJ_force(r) + elec_force(r, E_charge)
        r_array[i][z] = np.linalg.norm(r)
        z += 1 
        if j > i: 
            pres += np.linalg.norm(r) * gradU(np.linalg.norm(r), E_charge)
    return F, r_array, pres/2

def LJ_force(r):
    """F = - grad*U: the force resulting from the Lennard-Jones potential
    Args:
        r (array): the distance between each pair of particles
    Returns:
        F (array): the force between the particles
    """
    F = 4*(12*np.linalg.norm(r)**(-13) -6*np.linalg.norm(r)**(-7))*r/np.linalg.norm(r)
    return F

def elec_force(r, E_charge):
    """The Coulomb force between the particles
    Args:
        r (array): the distance between each pair of particles
        E_field (constant): the electrical field strength
    Returns:
        F (array): the Coulomb force"""
    F = -(k_c * E_charge**2 /(np.linalg.norm(r)**2))*r/np.linalg.norm(r) #Coulomb force
    F[2] += E_charge*E_field/N_particles
    return F

def gradU(r, E_charge):
  """ dU/dr
  Args:
      r (array): The distance between each pair of particles.
  Returns:
      dU (scalar): The gradient of the potential."""
  dU = 4*(-12*np.linalg.norm(r)**(-13) + 6*np.linalg.norm(r)**(-7))
  dU += k_c * E_charge**2 /(np.linalg.norm(r)**2)
  dU += -E_charge*E_field
  return dU

def potentialenergy(pos):
    """Calculates the total Lennard-Jones + Coulomb potential energy
    Args:
        pos (array): The position at the current timestep.
    Returns:
        U_LJ (scalar): The Lennard-Jones potential energy.
        U_C (scalar): The Coulomb potential energy."""
    U_LJ = 0 
    U_C = 0
    for i in range(N_particles):
      for j in range(i+1, N_particles):
        r = np.linalg.norm((pos[i] - pos[j] + L_sys/2)%L_sys - L_sys/2)
        U_LJ += potential(r, E_charge)[0]
        U_C += potential(r, E_charge)[1]
    return U_LJ, U_C

def elec_potentialenergy(E_charge, delta_r):
    """Calculates the electrical field contribution to the potential energy.
    Args:
        E-charge (scalar): the charge of the particles
        delta_r (array): the distance between the particles
    Returns:
        U (scalar): the electric field potential energy"""
    U = -np.sum(E_charge*E_field*delta_r[:,2])
    return U

def set_velocity():
    """Generates the initial velocities."""
    vel = np.random.normal(0, sigma, (N_particles, Dim))
    return vel

def scalefactor(kin):
    """Calculates the scale factor of the kinetic energy compared to the equilibrium value"""
    return np.sqrt(((3/2)*(N_particles-1)*Temp)/kin)

def kineticenergy(vel):
    """To calculate the total kinetic energy of the system.
    Args:
        vel (array): the velocities of the particles"""
    return np.sum(np.linalg.norm(vel, axis = 1)**2)/2 

def potential(r, E_charge):
    """The potential energy in a Lennard Jones potential
    Args:
        r (array): The distance between each pair of particles.
        E_charge (scalar): The electric charge of the particles.
    Returns:
        U_LJ (scalar): The potential energy from the Lennard Jones potential.
        U_C (scalar): The Coulombic potential energy."""
    U_LJ = 4*(r**(-12)-r**(-6)) 
    U_C = k_c * E_charge**2 /np.linalg.norm(r)
    return U_LJ, U_C

def generatelattice(): 
    """To generate the initial positions/velocities of the particles in the FCC lattice
    3x3x3 unit cells.
    Returns:
        pos (array): the initial positions."""
    a = L_sys/3
    pos = np.zeros((N_particles, Dim))
    c = 0
    for x in range(3):
      for y in range(3):
        for z in range(3):
          pos[4*c] = [0 + a*x, 0 + a*y ,0 + a*z] 
          pos[4*c+1] = [a/2 + a*x, a/2 + a*y , 0 + a*z] 
          pos[4*c+2] = [a/2 + a*x, 0 + a*y , a/2 + a*z] 
          pos[4*c+3] = [0 + a*x, a/2 + a*y , a/2 + a*z]
          c += 1 #counter
    return pos

def savedata(filename, hist, hist2, kin, LJ_pot, elec_pot, N_particles, rho, Temp, E_field, pres, h):
    """To save the data as txt files.
    Args:
        filename
        hist, hist2 (array): the two histograms of the distance between the particles.
        kin (array): the total kinetic energy for each time step.
        pot (array): the total potential energy for each time step.
        N_particles, rho, Temp, E_field, h (constants): the settings that were used for the simulation.
        pres (scalar): the pressure."""
    np.savetxt('histogram'+filename, (hist[0]+hist2[0])/2, delimiter =',')
    np.savetxt('bins'+filename, hist[1], delimiter = ',')
    np.savetxt('kineticenergy'+filename, kin, delimiter = ',')
    np.savetxt('LJ_potential'+filename, LJ_pot, delimiter = ',')
    np.savetxt('elec_potential'+filename, elec_pot, delimiter = ',')
    
    f = open('simulationparameters'+filename, 'x')
    f.write('Number of particles: {0} \n'.format(N_particles))
    f.write('Density {0} \n'.format(rho))
    f.write('Temperature: {0} \n'.format(Temp))
    f.write('E_field: {0} \n'.format(E_field))
    f.write('Pressure: {0} \n'.format(pres))
    f.write('Step size: {0}'.format(h))
    f.close()

def savetraj(data):
    """Source: https://stackoverflow.com/questions/3685265/how-to-write-a-multidimensional-array-to-a-text-file
    To save the trajectories of the particles."""
    with open('traj'+filename, 'w') as outfile:
        outfile.write('# Array shape: {0}\n'.format(data.shape))
        for data_slice in data:  
            np.savetxt(outfile, data_slice, fmt='%-7.3f')
            outfile.write('# New slice\n')
    pass


#The initial positions and velocities of the simulation:
pos_i = generatelattice()
vel_i = set_velocity()
#Empty arrays to save the data:
trajectories = np.zeros((N_particles, sim_time, Dim))
kin = np.zeros(sim_time)
LJ_pot = np.zeros(sim_time)
elec_pot = np.zeros(sim_time)
pres_array = np.zeros(sim_time)


l = 10 #to get the while loop started
while np.abs(l - 1) > threshold:
  for i in range(test_time):
    pos_i, vel_i, r, pres_array[i] = moveparticles(pos_i, vel_i, 0)
    trajectories[:,i,:] = pos_i
    delta_r = trajectories[:,i,:] - trajectories[:,i-1,:]
    kin[i] = kineticenergy(vel_i)
    LJ_pot[i], U_C = potentialenergy(pos_i)
    elec_pot[i] += U_C
    elec_pot[i:] += elec_potentialenergy(0, delta_r)
  l = scalefactor(kin[i])
  vel_i = vel_i*l
 
#now the actual simulation
trajectories = np.zeros((N_particles, sim_time, Dim))
kin = np.zeros(sim_time)
LJ_pot = np.zeros(sim_time)
elec_pot = np.zeros(sim_time)
pres_array = np.zeros(sim_time)
for i in range(sim_time):
  pos_i, vel_i, r, pres_array[i] = moveparticles(pos_i, vel_i, E_charge)
  trajectories[:,i,:] = pos_i
  delta_r = (trajectories[:,i,:] - trajectories[:,i-1,:] +L_sys/2)%L_sys - L_sys/2
  kin[i] = kineticenergy(vel_i)
  LJ_pot[i], U_C = potentialenergy(pos_i)
  elec_pot[i] += U_C
  elec_pot[i:] += elec_potentialenergy(E_charge, delta_r)
  if i == int(sim_time/2): 
      hist = np.histogram(r, bins)
  if i == sim_time-1:
      hist2 = np.histogram(r, bins)
   
pressure = (1 - ( np.mean(pres_array) / (3*N_particles*Temp) ))*(Temp*rho)
savedata(filename, hist, hist2, kin, LJ_pot, elec_pot, N_particles, rho, Temp, E_field, pressure, h)
savetraj(trajectories)
