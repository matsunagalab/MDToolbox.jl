#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_153.in \
 -c ../1_eq/run.rst \
 -o run_153.out \
 -r run_153.rst \
 -x run_153.nc

