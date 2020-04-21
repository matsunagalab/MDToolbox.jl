#!/bin/bash

program_at=/Users/yasu/genesis/bin/atdyn
export OMP_NUM_THREADS=1
name=C7eq
mpirun -np 8 ${program_at} ${name}.inp > ${name}.log
name=C7ax
mpirun -np 8 ${program_at} ${name}.inp > ${name}.log
