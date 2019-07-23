#!/bin/bash

NPROC=8
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_-105_105.in \
 -c ../3_eq/run_-105_105.rst \
 -o run_-105_105.out \
 -r run_-105_105.rst \
 -x run_-105_105.nc

