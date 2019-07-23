#!/usr/bin/python
# coding: utf-8

############# define functions
def print_lines(f, phi, psi):
    f.write("harmonic restraint\n")
    f.write(" &rst\n")
    f.write("   iat=5,7,9,15,\n")
    f.write("   r0=%f, k0=50.0,\n" % (phi))
    f.write(" /\n")
    f.write(" &rst\n")
    f.write("   iat=7,9,15,17,\n")
    f.write("   r0=%f, k0=50.0,\n" % (psi))
    f.write(" /\n")
    f.write("\n")

############# main
phi = range(-180, 0, 15)
psi = range(0, 180, 15)

for i in phi:
    for j in psi:
        filename = "run_%d_%d.disang" % (i, j)
        print "writing %s..." % (filename)
        f = open(filename, 'w')
        print_lines(f, i, j)
        f.close()

