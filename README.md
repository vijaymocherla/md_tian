# md_tian

`md_tian` (Molecular Dynamics Xia Tian) is a MD program written in Fortran for simulating the scattering of atoms (and molecules) from a surface. 

## Source code 

Following are a few details about the source code.

**List of modules:**

- `atom_class.f90` contains definitions of user types and all constants
- `force.f90` contains code to calculate energy and forces
- `mdalgo.f90` contains propagation algorithms
- `md_init.f90` contains code to set up a simulated system
- `md_tian.f90` is a main file governing simulations
- `open_file.f90` contains routines to open files smoothly
- `output.f90` contains output routines
- `useful_things.f90`	contains useful math routines

**Input files:**

- `md_tian.inp`	contains control parameters defining the simulation conditions
- `emt_{pes}.nml` contain emt-parameters for a species

## Installation

Compilation and linking (Intel Fortran):
```
$ ifort -O3 -ipo -o md_tian atom_class.f90 open_file.f90 useful_things.f90 md_init.f90 output.f90 mdalgo.f90 force.f90 md_tian.f90
```

The 1st working and tested version is put together February 18, 2014 
on a Fassberg Hill in Dynamics at Surfaces Dep. of MPIbpc
to the flaming storm of applause muffled by thick institute building walls.

## Credits:

>**Svenja Maria Janke**  
**Sascha Kandratsenka**  
**Daniel J. Auerbach**
>
>**Department of Dynamics at Surfaces**  
MPI for Biophysical Chemistry  
Am Fassberg 11  
37077 Goettingen  
Germany
>
>**Dynamics at Surfaces Department**  
Institute for Physical Chemistry  
Tammannstr. 6  
37077 Goettingen  
Germany

*`md_tian` is a very important program. It helps to better the world.*   

*jqrw sxrw=n! wr wj nA r sxrw nTrw!*
