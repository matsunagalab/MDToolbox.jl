#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_69.in \
 -c ../2_eq/run_69.rst \
 -o run_69.out \
 -r run_69.rst \
 -x run_69.nc

