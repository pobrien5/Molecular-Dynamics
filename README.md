# Molecular Dynamics of a Argon atoms
-------------------------------------

Date: 24/03/22
Leiden University

Authors (Group 14): Marieke Visscher \& Patrick O'Brien

-------------------------------------

Thank you for downloading our code.
Instructions on how to use it are given below.

(0) Unzip the compressed file into a location of your choice.
This file contains two python scripts - One to run the simulation
(MD_simulation.py) and one to plot the data (MD_plots.py).
(1) Launch MD_simulation.py with your favourite python environment.
(2) Set your simulation parameters...

	line 12: Temp = ?? (float)
	line 13: rho = ?? (float)
	line 14: filename = ?? + '.txt' (string)

Temp (temperature) and rho (number density) determine the
phase of your system. Both quantities are in dimensionless
units. The filename determines what your data will be 
called when it is saved - Entirely up to you!
Other parameters you might want to change are...

	line 16: sim_time = ?? (integer)
	line 17: E_field = ?? (float)
	line 18: E_charge = ?? (integer)

By default these are 100, 0.0 and 0 respectively. Try setting 
E_field = 5.0 and E_charge = 1 for example. 
(3) Run the code. This will take ~10 minutes for 100 timesteps.
The code saves the data in 7 text files... 

* 'traj'+filename+'.txt' - positions (x,y,z) of each particle at
   each timestep.
* 'simulationparameters'+filename+'.txt' - saves important
   information about this specific simulation. 
* 'bins'+filename+'.txt' - includes the bins of the histogram
* 'histogram'+filename+'.txt' - includes the histogram of the 
   distance between particles. 
* 'kineticenergy'+filename+'.txt' - kinetic energy of the sytem
   over the course of the simulation.
* 'LJ_potential'+filename+'.txt' - Lennard-Jones potential energy
   of the sytem over the course of the simulation.
* 'elec_potential'+filename+'.txt' - Coulombic + electric field
   potential energy of the simulation.

(4) Launch MD_plots.py to analyse this data.
(5) Fill in your simulation parameters from before...

	line 15: rho = ?? (float)
	line 16: sim_time = ?? (integer)

(6) Run MD_plots.py to generate the plots + play the movie.







