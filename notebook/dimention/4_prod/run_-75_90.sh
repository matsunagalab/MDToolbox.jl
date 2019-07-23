#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-75_90.in \
 -c ../3_eq/run_-75_90.rst \
 -o run_-75_90.out \
 -r run_-75_90.rst \
 -x run_-75_90.nc

