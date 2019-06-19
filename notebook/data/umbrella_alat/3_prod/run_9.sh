#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_9.in \
 -c ../2_eq/run_9.rst \
 -o run_9.out \
 -r run_9.rst \
 -x run_9.nc

