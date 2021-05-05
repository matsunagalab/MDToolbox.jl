#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_165.in \
 -c ../1_eq/run.rst \
 -o run_165.out \
 -r run_165.rst \
 -x run_165.nc

