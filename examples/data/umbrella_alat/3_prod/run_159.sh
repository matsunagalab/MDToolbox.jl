#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_159.in \
 -c ../2_eq/run_159.rst \
 -o run_159.out \
 -r run_159.rst \
 -x run_159.nc

