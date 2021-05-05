#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_78.in \
 -c ../1_eq/run.rst \
 -o run_78.out \
 -r run_78.rst \
 -x run_78.nc

