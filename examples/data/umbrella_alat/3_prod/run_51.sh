#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_51.in \
 -c ../2_eq/run_51.rst \
 -o run_51.out \
 -r run_51.rst \
 -x run_51.nc

