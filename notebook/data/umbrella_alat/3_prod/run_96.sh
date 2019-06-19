#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_96.in \
 -c ../2_eq/run_96.rst \
 -o run_96.out \
 -r run_96.rst \
 -x run_96.nc

