#!/usr/bin/python
# coding: utf-8

import os

############# define functions
def print_lines(f, phi, psi):
    f.write("#!/bin/bash\n")
    f.write("\n")
    f.write("NPROC=8\n")
    f.write("mpirun -np $NPROC sander.MPI -O \\\n")
    #f.write("pmemd.cuda -O \\\n")
    f.write(" -p prmtop \\\n")
    f.write(" -i run_%d_%d.in \\\n" % (phi, psi))
    f.write(" -c ../3_eq/run_%d_%d.rst \\\n" % (phi, psi))
    f.write(" -o run_%d_%d.out \\\n" % (phi, psi))
    f.write(" -r run_%d_%d.rst \\\n" % (phi, psi))
    f.write(" -x run_%d_%d.nc\n" % (phi, psi))
    f.write("\n")

############# main
phi = range(-180, 0, 15)
psi = range(0, 180, 15)

for i in phi:
    for j in psi:
        filename = "run_%d_%d.sh" % (i, j)
        print "writing %s..." % (filename)
        f = open(filename, 'w')
        print_lines(f, i, j)
        f.close()
        os.chmod(filename, 0755)

