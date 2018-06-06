#!/usr/bin/python

# Script for generating plots using matplotlib in python
# Anything after a '#' character is ignored

import argparse
import numpy as np
import os
import subprocess
import math
import sys

import pvparse as pvp

fig_num = 0
PROGRESS = 10
WIDTH = 80

pargs = ['b', 'r']

parser = argparse.ArgumentParser()

parser.add_argument("file", help="name of file containing position information")
parser.add_argument("-o", "--out", help="root name of output graphics/movie")
parser.add_argument("-x", "--xlim", help="x-axis limits on plots",
                    nargs='+', type=float, default=1.0)
parser.add_argument("-y", "--ylim", help="y-axis limits on plots",
                    nargs='+', type=float, default=1.0)
parser.add_argument("-p", "--pad", type=int, default=3,
                    help="number of digits in last plot")
parser.add_argument("-m", "--movie", help="build an mp4 instead of the plots",
                    action="store_true")
parser.add_argument("-w", "--warmstart", help="skips the image building step",
                    action="store_true")
parser.add_argument("-k", "--keep", help="keep position plot files",
                    action="store_true")

args = parser.parse_args()

if (args.out == None):
    root, sep, last = str(args.file).rpartition(".")
else:
    root = args.out
format_str = "_%0" + str(args.pad) + "d.png"

if(not args.warmstart):
    print "Parsing file and building plots..."

    import matplotlib.pyplot as plt

    if(isinstance(args.xlim, float)):
        a = -abs(args.xlim)
        b = abs(args.xlim)
    elif(len(args.xlim) < 2):
        a = -abs(args.xlim[0])
        b = abs(args.xlim[0])
    else:
        a = args.xlim[0]
        b = args.xlim[1]
        if (a >= b):
            print "Error: Left limit", a, "must be less than right limit", b
            quit()

    if(isinstance(args.ylim, float)):
        c = -abs(args.ylim)
        d = abs(args.ylim)
    elif(len(args.ylim) < 2):
        c = -abs(args.ylim[0])
        d = abs(args.ylim[0])
    else:
        c = args.ylim[0]
        d = args.ylim[1]
        if (c >= d):
            print "Error: Lower limit", c, "must be less than upper limit", d
            quit()

    data = open(args.file, mode='r')

    nentries = pvp.findNumEntries(data)
    if (nentries < 1):
        print "Number of entries is non positive -- got", nentries
        quit()
    else:
        pad = len(str(nentries))
        format_str = "_%0" + str(pad) + "d.png"

    processed = 0
    while (processed < nentries):

        pvp.reportProgress(processed, nentries, WIDTH)

        sts, num, dim, time, layer, vort, pos = pvp.readPositionEntry(data)

        if (sts != 0):
            print "**ERROR:: readPositionEntry returned", sts
            quit()

        # Create and save the plot for this entry --
        # this won't work on sphere
        lys = np.unique(layer)
        if lys.size == 1:
            fig = plt.figure(figsize=(8,8))
            plt.scatter(pos[:,0], pos[:,1], c=pargs[0])
            plt.xlim(a, b)
            plt.ylim(c, d)
            plt.legend()
        else:
            fig, axes = plt.subplots(nrows=1, ncols=lys.size,
                                     figsize=(8*lys.size, 8))
            for k in range(lys.size):
                indx = np.nonzero(layer == lys[k])
                # axes[k].scatter(pos[indx,0], pos[indx,1], c=pargs[k],
                #             label='Layer ' + str(pvp.LAYER_MAP[lys[k]]))
                # axes[k].legend()
                axes[k].scatter(pos[indx,0], pos[indx,1], c=pargs[k])
                axes[k].set_title('Layer ' + str(pvp.LAYER_MAP[lys[k]]))
                axes[k].set_xlim(a, b)
                axes[k].set_ylim(c, d)

        title = "Simulation Time:%8.3f" % time
        plt.suptitle(title)
        file_str = root + format_str % fig_num
        plt.savefig(file_str)
        plt.close()
        fig_num = fig_num + 1

        processed += 1

    data.close()
    if (processed != nentries):
        print "WARN: Number of expected entries", nentries
        print "Found", processed
    else:
        pvp.reportProgress(nentries, nentries, WIDTH)

if(args.movie):
    print "Building movie..."
    files = root + format_str
    movie = root + ".mp4"
    if os.path.exists(movie):
        os.remove(movie)
    pargs = ["ffmpeg", "-i", files, "-c:v", "libx264", "-pix_fmt",
             "yuv420p", movie]
    subprocess.call(pargs)
    if (not args.warmstart and not args.keep):
        # Remove image files
        for k in range(fig_num):
            file_str = root + format_str % k
            os.remove(file_str)
