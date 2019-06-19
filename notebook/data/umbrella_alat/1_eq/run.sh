#!/bin/bash

NPROC=2
mpirun -np $NPROC sander.MPI -O \
 -p prmtop \
 -i run.in \
 -c crd \
 -o run.out \
 -r run.rst \
 -x run.nc

