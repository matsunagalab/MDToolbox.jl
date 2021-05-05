#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_21.in \
 -c ../1_eq/run.rst \
 -o run_21.out \
 -r run_21.rst \
 -x run_21.nc

