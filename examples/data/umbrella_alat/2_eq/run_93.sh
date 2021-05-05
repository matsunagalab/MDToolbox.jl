#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_93.in \
 -c ../1_eq/run.rst \
 -o run_93.out \
 -r run_93.rst \
 -x run_93.nc

