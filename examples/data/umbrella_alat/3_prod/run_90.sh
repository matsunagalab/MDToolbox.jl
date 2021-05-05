#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_90.in \
 -c ../2_eq/run_90.rst \
 -o run_90.out \
 -r run_90.rst \
 -x run_90.nc

