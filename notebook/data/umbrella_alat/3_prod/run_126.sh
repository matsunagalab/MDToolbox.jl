#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run_126.in \
 -c ../2_eq/run_126.rst \
 -o run_126.out \
 -r run_126.rst \
 -x run_126.nc

