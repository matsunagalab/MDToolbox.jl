#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_156.in \
 -c ../2_eq/run_156.rst \
 -o run_156.out \
 -r run_156.rst \
 -x run_156.nc

