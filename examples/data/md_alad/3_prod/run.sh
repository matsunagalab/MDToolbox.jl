#!/bin/bash -x

pmemd.cuda -O \
 -p ../0_init/alad_water.prmtop \
 -i run.in \
 -c ../2_eq/run.rst \
 -o run.out \
 -r run.rst \
 -x run.nc


