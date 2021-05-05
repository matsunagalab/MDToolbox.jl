#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_15.in \
 -c ../1_eq/run.rst \
 -o run_15.out \
 -r run_15.rst \
 -x run_15.nc

