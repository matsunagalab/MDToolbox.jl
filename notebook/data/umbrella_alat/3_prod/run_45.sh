#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_45.in \
 -c ../2_eq/run_45.rst \
 -o run_45.out \
 -r run_45.rst \
 -x run_45.nc

