#!/usr/bin/python

# Script for building plots of observables using optimal closure theory

import argparse
import numpy as np
import scipy.special

import pvparse as pvp
import poisson_solver as ps

import os

GRID_SIZE = 128
WIDTH = 80

# Set locations for various things
home = os.path.expanduser("~")
OUTPUT = home + "/Vortex/output/"
MOVIE = home + "/Vortex/movies/"

parser = argparse.ArgumentParser()

parser.add_argument("root", help="root of output file name to get results from")
parser.add_argument("-n", "--nentries", help="number of entries per file",
                    required=True, type=int)
parser.add_argument("-r", "--nfiles", help="number of files in batch",
                    required=True, type=int)
parser.add_argument("-v", "--npvs",
                    help="number of point vortices in simulation",
                    type=int, required=True)
parser.add_argument("-d", "--defrad", type=float, required=True,
                    help="deformation radius")

args = parser.parse_args()

nentries = args.nentries
nfiles = args.nfiles
npvs = args.npvs
rd = args.defrad
kd = np.sqrt(2.0) / rd
gamma = 2.0 / npvs

posfmt = MOVIE + args.root + "/" + args.root + "_%0"
posfmt += str(len(str(nfiles))) + "d.pos"

#grid = ps.Grid(GRID_SIZE)
grid = ps.LegGrid(GRID_SIZE, 4*rd)
APEavg = np.zeros(nentries)

for l in range(nfiles):
    posfile = posfmt % l
    print "Processing file", posfile
    print "Reading file..."
    sts, ti, layer, vort, pos = pvp.readPositionFile(posfile, npvs, 2, True)
    print "Computing..."
    for i in range(nentries):
        pvp.reportProgress(i, nentries, WIDTH)
        dx = (pos[i,:,0] -
              np.tile(grid.xmat, (npvs,1,1)).transpose()).transpose()
        dy = (pos[i,:,1] -
              np.tile(grid.ymat, (npvs,1,1)).transpose()).transpose()
        # ri = kd * np.sqrt(dx*dx + dy*dy)
        ri = kd * np.sqrt(dx*dx + dy*dy)
        k0 = scipy.special.k0(ri)
        psic = (-layer * k0.transpose()).transpose()
        # psic = (-layer * k0.transpose()).transpose()
        psic = gamma * 0.25 / np.pi * np.sum(psic, 0)
        ape = kd**2 * np.sum(grid.wmat * psic**2)
        APEavg[i] += 1.0 / (l + 1.) * (ape - APEavg[i])

    pvp.reportProgress(nentries, nentries, WIDTH)

outfile = OUTPUT + args.root + "/" + args.root + ".ape"
out = open(outfile, 'w')

for k in range(APEavg.size):
    output = "time={:.15e}".format(ti[k])
    output += "  ape={:.15e}".format(APEavg[k])
    output += "\n"
    out.write(output)

out.flush()
out.close()
