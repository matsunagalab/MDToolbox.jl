#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_135.in \
 -c ../1_eq/run.rst \
 -o run_135.out \
 -r run_135.rst \
 -x run_135.nc

