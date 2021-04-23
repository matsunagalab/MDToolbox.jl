#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_132.in \
 -c ../2_eq/run_132.rst \
 -o run_132.out \
 -r run_132.rst \
 -x run_132.nc

