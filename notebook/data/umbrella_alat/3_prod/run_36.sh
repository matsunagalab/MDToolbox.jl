#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_36.in \
 -c ../2_eq/run_36.rst \
 -o run_36.out \
 -r run_36.rst \
 -x run_36.nc

