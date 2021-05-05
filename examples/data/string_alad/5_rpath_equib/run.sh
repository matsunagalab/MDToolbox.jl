#!/bin/bash

program_at=/Users/yasu/genesis/bin/atdyn
export OMP_NUM_THREADS=1
mpirun -np 16 ${program_at} rpath.1.inp | tee rpath.1.log

