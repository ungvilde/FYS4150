# Simulating particles in a Penning trap

We simulate the movement of particles in a Penning trap. 
The repository contains the following files and folders:

- **simulate_one_two_particles.cpp**: Run experiments for single particle and two particles in Penning trap
- **simulate_time_dependence**: Run experiments with time-varying electric field
- **plot.py**: plots the relevant figures
- **plot_functions.py**: holds the functions used for plotting
- **data/**: where data files from numerical experiments are stored
- **figs/**: where we have the plots

To compile and run the C++ programs, do

> `` make all``

Note that the simulations will take some time to complete, approximately 1-1.5h for all the numerical experiments.

## Particle and PenningTrap classes

The simulation is based on two classes: ``Particle`` and ``PenningTrap``. A ``Particle`` object holds the mass, charge, velocity and position of the particle, as well as methods for updating the position and velocity of the particle. A ``PenningTrap`` has a given configuration for the magnetic and electric fields, and it holds particles. Particles can be added and counted. There are also methods for updating the velocities and positions of the particles with one time step. The methods implemented are 4th-order Runge-Kutta and the Euler-Cromer method.
