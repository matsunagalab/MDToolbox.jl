#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_132.in \
 -c ../1_eq/run.rst \
 -o run_132.out \
 -r run_132.rst \
 -x run_132.nc

