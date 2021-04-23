#!/bin/bash

program_at=/Users/yasu/genesis/bin/atdyn
export OMP_NUM_THREADS=1
mpirun -np 16 ${program_at} rpath.2.inp | tee rpath.2.log

