#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_180.in \
 -c ../2_eq/run_180.rst \
 -o run_180.out \
 -r run_180.rst \
 -x run_180.nc

