#!/bin/bash

rm *.mod
#ifort -o md_tian_0003 -mkl=sequential open_file.f90 atom_class.f90 useful_things.f90 md_init.f90 force.f90 dtrnlspbc.f90 fit4_tian.f90 mdalgo.f90 output.f90 md_tian.f90
#ifort -o md_tian_160701 -mkl=sequential open_file.f90 atom_class.f90 useful_things.f90 md_init.f90 force.f90 fit4_tian.f90 mdalgo.f90 output.f90 md_tian.f90
ifx -o md_tian.x -qmkl open_file.f90 atom_class.f90 useful_things.f90 md_init.f90 force.f90 dtrnlspbc.f90 fit4_tian.f90 mdalgo.f90 output.f90 md_tian.f90 
