#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_42.in \
 -c ../1_eq/run.rst \
 -o run_42.out \
 -r run_42.rst \
 -x run_42.nc

