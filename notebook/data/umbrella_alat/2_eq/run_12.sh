#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_12.in \
 -c ../1_eq/run.rst \
 -o run_12.out \
 -r run_12.rst \
 -x run_12.nc

