#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-150_135.in \
 -c ../3_eq/run_-150_135.rst \
 -o run_-150_135.out \
 -r run_-150_135.rst \
 -x run_-150_135.nc

