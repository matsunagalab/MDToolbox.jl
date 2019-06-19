#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_63.in \
 -c ../2_eq/run_63.rst \
 -o run_63.out \
 -r run_63.rst \
 -x run_63.nc

