#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-120_75.in \
 -c ../3_eq/run_-120_75.rst \
 -o run_-120_75.out \
 -r run_-120_75.rst \
 -x run_-120_75.nc

