#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_84.in \
 -c ../2_eq/run_84.rst \
 -o run_84.out \
 -r run_84.rst \
 -x run_84.nc

