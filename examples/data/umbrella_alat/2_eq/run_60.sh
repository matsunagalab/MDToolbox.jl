#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_60.in \
 -c ../1_eq/run.rst \
 -o run_60.out \
 -r run_60.rst \
 -x run_60.nc

