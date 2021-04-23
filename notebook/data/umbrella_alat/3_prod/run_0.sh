#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_0.in \
 -c ../2_eq/run_0.rst \
 -o run_0.out \
 -r run_0.rst \
 -x run_0.nc

