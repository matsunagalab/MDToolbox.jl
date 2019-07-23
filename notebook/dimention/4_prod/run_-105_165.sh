#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-105_165.in \
 -c ../3_eq/run_-105_165.rst \
 -o run_-105_165.out \
 -r run_-105_165.rst \
 -x run_-105_165.nc

