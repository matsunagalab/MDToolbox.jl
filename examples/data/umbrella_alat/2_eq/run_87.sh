#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_87.in \
 -c ../1_eq/run.rst \
 -o run_87.out \
 -r run_87.rst \
 -x run_87.nc

