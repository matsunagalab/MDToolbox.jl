#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_27.in \
 -c ../1_eq/run.rst \
 -o run_27.out \
 -r run_27.rst \
 -x run_27.nc

