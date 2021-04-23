#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_153.in \
 -c ../2_eq/run_153.rst \
 -o run_153.out \
 -r run_153.rst \
 -x run_153.nc

