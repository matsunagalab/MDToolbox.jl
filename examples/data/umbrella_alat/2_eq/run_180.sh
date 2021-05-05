#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_180.in \
 -c ../1_eq/run.rst \
 -o run_180.out \
 -r run_180.rst \
 -x run_180.nc

