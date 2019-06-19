#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_6.in \
 -c ../2_eq/run_6.rst \
 -o run_6.out \
 -r run_6.rst \
 -x run_6.nc

