#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_3.in \
 -c ../2_eq/run_3.rst \
 -o run_3.out \
 -r run_3.rst \
 -x run_3.nc

