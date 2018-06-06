#!/usr/bin/python

import argparse
import glob
import os
import sys

parser = argparse.ArgumentParser()

parser.add_argument("root", help="root of output file name to get results from")
parser.add_argument("-o", "--old-size", help="number of files in batch",
                    required=True, type=int)
parser.add_argument("-n", "--new-size", required=True, type=int,
                    help="new number of files in batch")
parser.add_argument("-d", "--directories", default=None, nargs="+",
                    help="directories in which things need to be renumbered",
                    required=True)

args = parser.parse_args()

root = args.root
nnum = args.new_size
onum = args.old_size
dirs = args.directories

onstr = "{:0" + str(len(str(onum))) + "d}"
nnstr = "{:0" + str(len(str(nnum))) + "d}"

if onum > nnum:
    print "Current batch size is bigger than new batch size.  Aborting."
    sys.exit(-1)

for cdir in dirs:
    for k in range(onum):
        fnum = onstr.format(k)
        nfnum = nnstr.format(k)
        fpattern = cdir + "/" + root + "/" + root + "_" + fnum + ".*"
        torename = glob.glob(fpattern)
        for cfile in torename:
            beg, n, rem = cfile.rpartition(fnum)
            fname = beg + nfnum + rem
            os.rename(cfile, fname)
