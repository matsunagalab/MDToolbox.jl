#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_141.in \
 -c ../1_eq/run.rst \
 -o run_141.out \
 -r run_141.rst \
 -x run_141.nc

