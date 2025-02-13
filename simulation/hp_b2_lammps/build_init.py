import math, argparse, random
import numpy as np

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('ntypes', type=int, help="number of types")
    parser.add_argument('s1', type=str, help="first chain sequence")
    parser.add_argument('s2', type=str, help="second chain sequence")
    parser.add_argument('--Lscale', type=float, default=4., help="scaling factor for box dimension [4]")
    parser.add_argument('filename', type=str, help='filename')
    clargs = parser.parse_args()

    filename = clargs.filename

    natoms = len(clargs.s1) + len(clargs.s2)
    nbonds = natoms - 2

    X0 = natoms * clargs.Lscale / 2.

    print("Writing to polymer data file...")
    with open(filename, 'w') as f:
        f.write("Sequences: %s,%s\n\n" % (clargs.s1, clargs.s2))
        f.write("%16d atoms\n" % natoms)
        f.write("%16d bonds\n" % nbonds)
        f.write("\n%16d atom types\n%16d bond types\n\n" % (clargs.ntypes, 1))
        for t in 'xyz':
            f.write("%16f %16f %slo %shi\n" % (0, 2. * X0, t, t))
        f.write("\nMasses\n\n")
        for i in range(clargs.ntypes):
            f.write("%16d %16d\n" % (i+1, 1.))
        f.write("\nAtoms # bond\n\n")
        for i in range(len(clargs.s1)):
            f.write("%16d %16d %16d %16f %16f %16f\n" % \
                    (i + 1, 1, int(clargs.s1[i]), X0 - len(clargs.s1) + i, X0, X0))
        for i in range(len(clargs.s2)):
            f.write("%16d %16d %16d %16f %16f %16f\n" % \
                    (len(clargs.s1) + i + 1, 2, int(clargs.s2[i]), X0 + i, X0, X0))
        f.write("\nBonds\n\n")
        for i in range(len(clargs.s1) - 1):
            f.write("%16d %16d %16d %16d\n" % (i + 1, 1, i + 1, i + 2))
        for i in range(len(clargs.s1) - 1, natoms - 2):
            f.write("%16d %16d %16d %16d\n" % (i + 1, 1, i + 2, i + 3))
