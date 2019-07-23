#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-105_60.in \
 -c ../3_eq/run_-105_60.rst \
 -o run_-105_60.out \
 -r run_-105_60.rst \
 -x run_-105_60.nc

