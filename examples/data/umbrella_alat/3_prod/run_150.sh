#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_150.in \
 -c ../2_eq/run_150.rst \
 -o run_150.out \
 -r run_150.rst \
 -x run_150.nc

