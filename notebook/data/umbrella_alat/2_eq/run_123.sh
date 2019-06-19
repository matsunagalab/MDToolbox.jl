#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_123.in \
 -c ../1_eq/run.rst \
 -o run_123.out \
 -r run_123.rst \
 -x run_123.nc

