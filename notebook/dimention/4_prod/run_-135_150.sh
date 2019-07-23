#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-135_150.in \
 -c ../3_eq/run_-135_150.rst \
 -o run_-135_150.out \
 -r run_-135_150.rst \
 -x run_-135_150.nc

