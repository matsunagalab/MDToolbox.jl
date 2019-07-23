#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-120_0.in \
 -c ../3_eq/run_-120_0.rst \
 -o run_-120_0.out \
 -r run_-120_0.rst \
 -x run_-120_0.nc

