#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-15_15.in \
 -c ../3_eq/run_-15_15.rst \
 -o run_-15_15.out \
 -r run_-15_15.rst \
 -x run_-15_15.nc

