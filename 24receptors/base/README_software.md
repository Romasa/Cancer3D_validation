CANCER3D
========

CANCER3D is a modelling framework devoted to simulating 3D tumor growth across several scales. It captures the microphysiology of individual cells and describes tissue mechanics with an agent-based approach. Brownian dynamics (BD) are used to preserve the stochastic and spatial aspects characterizing the regulation of single-cells.


# Dependencies

The CANCER3D software was developed and tested under Ubuntu 16.04. There exist few dependencies that need to be installed on the OS before compiling and running CANCER3D:

- GCC: http://gcc.gnu.org/
	- apt-get install gcc

- GNU Scientific Library (GSL) https://www.gnu.org/software/gsl/
	- download and install the GSL package: ftp://ftp.gnu.org/gnu/gsl/

# Compiling and running

GNU make is used to facilitate the compiling and running of the code. In Makefile, replace "USERNAME" with your linux username in the paths. Then use the following commands:

- To compile: make all
- To execute: make run

# Post-processing and visualization

VTK files are generated in the particle_vtk folder upon running the program. We recommand using ParaView to properly visualize the simulations.
