#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_33.in \
 -c ../1_eq/run.rst \
 -o run_33.out \
 -r run_33.rst \
 -x run_33.nc

