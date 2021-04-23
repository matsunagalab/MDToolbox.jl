#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_147.in \
 -c ../1_eq/run.rst \
 -o run_147.out \
 -r run_147.rst \
 -x run_147.nc

