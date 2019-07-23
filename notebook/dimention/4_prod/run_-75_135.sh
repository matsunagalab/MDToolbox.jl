#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-75_135.in \
 -c ../3_eq/run_-75_135.rst \
 -o run_-75_135.out \
 -r run_-75_135.rst \
 -x run_-75_135.nc

