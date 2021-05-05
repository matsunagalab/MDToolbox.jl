#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_18.in \
 -c ../2_eq/run_18.rst \
 -o run_18.out \
 -r run_18.rst \
 -x run_18.nc

