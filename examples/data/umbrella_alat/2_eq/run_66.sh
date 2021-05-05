#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_66.in \
 -c ../1_eq/run.rst \
 -o run_66.out \
 -r run_66.rst \
 -x run_66.nc

