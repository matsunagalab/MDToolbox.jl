#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_48.in \
 -c ../1_eq/run.rst \
 -o run_48.out \
 -r run_48.rst \
 -x run_48.nc

