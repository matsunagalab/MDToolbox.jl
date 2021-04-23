#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_168.in \
 -c ../1_eq/run.rst \
 -o run_168.out \
 -r run_168.rst \
 -x run_168.nc

