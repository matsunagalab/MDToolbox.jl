#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-15_45.in \
 -c ../3_eq/run_-15_45.rst \
 -o run_-15_45.out \
 -r run_-15_45.rst \
 -x run_-15_45.nc

