#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_6.in \
 -c ../1_eq/run.rst \
 -o run_6.out \
 -r run_6.rst \
 -x run_6.nc

