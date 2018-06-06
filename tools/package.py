#!/usr/bin/python

import argparse
import os
import subprocess

# Set locations for various things
home = os.path.expanduser("~")
VORTEX = home + "/Vortex/"
SCRIPTS = home + "/Vortex/tools/"
RESULTS = home + "/Vortex/output/"
CLOSE = home + "/Vortex/results/"
LOG = home + "/Vortex/log/"
DATA = home + "/Vortex/data/"
BIN = home + "/Vortex/bin/"
MOVIE = home + "/Vortex/movies/"

# Useful lists
dirs = [DATA, RESULTS, MOVIE, CLOSE]
app = ["_in", "_out", "_pos", "_oc"]

# Create command line argument parser
parser = argparse.ArgumentParser()

# Define accepted arguments
parser.add_argument("root", help="root name for files produced by these runs")
parser.add_argument("-u", "--unpack", help="unpack the given archives",
                    action="store_true")
parser.add_argument("-s", "--skip-position", action="store_true",
                    help="skip the position files in packing or unpacking")

# Parse arguments
args = parser.parse_args()
bname = args.root

if args.unpack:
    targs = ["tar", "-xf"]
    gargs = ["gzip", "-d"]
    for i in range(len(dirs)):
        os.chdir(dirs[i])
        arc = VORTEX + bname + app[i]
        target = dirs[i] + bname
        if os.path.exists(arc + ".tar.gz"):
            if os.path.exists(target):
                print target, "already exists"
            else:
                if dirs[i] == MOVIE and args.skip_position:
                    print "Moving", arc + ".tar.gz", "to", MOVIE
                    os.rename(arc + ".tar.gz", MOVIE + arc + ".tar.gz")
                else:
                    print "Unpacking", arc + ".tar.gz"
                    pargs = gargs + [arc + ".tar.gz"]
                    retcode = subprocess.call(pargs)
                    pargs = targs + [arc + ".tar"]
                    retcode = subprocess.call(pargs)
                    os.remove(arc + ".tar")
        else:
            print "No archive for", target
else:
    targs = ["tar", "-cf"]
    gargs = ["gzip", "--best"]

    for i in range(len(dirs)):
        os.chdir(dirs[i])
        aname = bname + app[i] + ".tar"
        tname = bname
        if os.path.exists(tname):
            if dirs[i] == MOVIE and args.skip_position:
                print "Skipping", dirs[i] + bname + app[i] + ".tar.gz"
                continue
            else:
                print "Archiving", dirs[i] + tname
                pargs = targs + [aname, tname]
                retcode = subprocess.call(pargs)
                pargs = gargs + [aname]
                retcode = subprocess.call(pargs)
                os.rename(dirs[i] + bname + app[i] + ".tar.gz",
                          VORTEX + bname + app[i] + ".tar.gz")
        else:
            print "Not found in", dirs[i]
