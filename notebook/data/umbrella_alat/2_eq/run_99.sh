#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_99.in \
 -c ../1_eq/run.rst \
 -o run_99.out \
 -r run_99.rst \
 -x run_99.nc

