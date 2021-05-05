#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_171.in \
 -c ../1_eq/run.rst \
 -o run_171.out \
 -r run_171.rst \
 -x run_171.nc

