#!/usr/bin/python

import argparse
import numpy as np
import matplotlib.pyplot as plt

import pvparse as pvp

# Create command line argument parser
parser = argparse.ArgumentParser()

# Define accepted arguments
parser.add_argument("file", help="name of the file to write initial data to")
parser.add_argument("-n", "--number", help="number of vortices in the file",
                    type=int)
parser.add_argument("-d", "--dimension", help="dimension vortices live in",
                    type=int)

# Parse the command line
args = parser.parse_args()

num = args.number
dim = args.dimension

sts, vort, pos = pvp.readInitFile(args.file, num, dim)
if (sts != 0):
    print "Parsing file returned", sts
    quit()

if (dim == 2):
    plt.figure()
    plt.plot(pos[:,0], pos[:,1], 'bo')
    root, sep, end = args.file.rpartition(".")
    fstr = root + ".png"
    plt.savefig(fstr)
    plt.close()
    print "Saved plot as", fstr
else:
    print "Dimension", dim, "not implemented"
