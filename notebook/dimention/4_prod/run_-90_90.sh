#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-90_90.in \
 -c ../3_eq/run_-90_90.rst \
 -o run_-90_90.out \
 -r run_-90_90.rst \
 -x run_-90_90.nc

