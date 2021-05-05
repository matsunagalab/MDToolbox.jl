#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_57.in \
 -c ../1_eq/run.rst \
 -o run_57.out \
 -r run_57.rst \
 -x run_57.nc

