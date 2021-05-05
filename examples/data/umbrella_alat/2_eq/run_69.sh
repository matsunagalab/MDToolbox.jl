#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_69.in \
 -c ../1_eq/run.rst \
 -o run_69.out \
 -r run_69.rst \
 -x run_69.nc

