#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-165_165.in \
 -c ../3_eq/run_-165_165.rst \
 -o run_-165_165.out \
 -r run_-165_165.rst \
 -x run_-165_165.nc

