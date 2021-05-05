#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_33.in \
 -c ../2_eq/run_33.rst \
 -o run_33.out \
 -r run_33.rst \
 -x run_33.nc

