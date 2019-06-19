#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_48.in \
 -c ../2_eq/run_48.rst \
 -o run_48.out \
 -r run_48.rst \
 -x run_48.nc

