#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_174.in \
 -c ../1_eq/run.rst \
 -o run_174.out \
 -r run_174.rst \
 -x run_174.nc

