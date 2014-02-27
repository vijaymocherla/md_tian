md_tian
=======

md_tian (Molecular Dynamics Xia Tian) is a program for sumulating the scattering of atoms (and molecules) 
from a surface. 

Source code is in Fortran.
List of modules:

atom_class.f90		contains definitions of user types and all constants
force.f90		contains code to calculate energy and forces
mdalgo.f90		contains propagation algorithms
md_init.f90		contains code to set up a simulated system
md_tian.f90		is a main file governing simulations
open_file.f90		contains routines to open files smoothly
output.f90		contains output routines
useful_things.f90	contains useful math routines

Input files:

md_tian.inp	contains control parameters defining the simulation conditions
*.nml		contain emt-parameters for a species

Compilation and linking (Intel Fortran):

ifort -O3 -ipo -o md_tian atom_class.f90 open_file.f90 useful_things.f90 md_init.f90 output.f90 mdalgo.f90 force.f90 md_tian.f90


The 1st working and tested version is put together February 18, 2014 
on a Fassberg Hill in Dynamics at Surfaces Dep. of MPIbpc
to the flaming storm of applause muffled by thick institute building walls.

Credits:

Svenja Maria Janke	
Sascha Kandratsenka
Daniel J. Auerbach

Dynamics at Surfaces Dep.
MPI for biophysical Chemistry
Am Fassberg 11
37077 Goettingen
Germany

Dynamics at Surfaces Dep.
Institute for Physical Chemistry
Tammannstr. 6
37077 Goettingen
Germany

Md xia4 tian1 is a very important program. It helps to better the world.
jqrw sHrw=n! wr wj nA r sHrw nTrw!
