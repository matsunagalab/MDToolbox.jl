#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_18.in \
 -c ../1_eq/run.rst \
 -o run_18.out \
 -r run_18.rst \
 -x run_18.nc

