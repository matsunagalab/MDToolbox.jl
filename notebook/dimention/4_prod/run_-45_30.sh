#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-45_30.in \
 -c ../3_eq/run_-45_30.rst \
 -o run_-45_30.out \
 -r run_-45_30.rst \
 -x run_-45_30.nc

