#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_63.in \
 -c ../1_eq/run.rst \
 -o run_63.out \
 -r run_63.rst \
 -x run_63.nc

