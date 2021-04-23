#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_3.in \
 -c ../1_eq/run.rst \
 -o run_3.out \
 -r run_3.rst \
 -x run_3.nc

