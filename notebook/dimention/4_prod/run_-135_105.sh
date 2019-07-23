#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-135_105.in \
 -c ../3_eq/run_-135_105.rst \
 -o run_-135_105.out \
 -r run_-135_105.rst \
 -x run_-135_105.nc

