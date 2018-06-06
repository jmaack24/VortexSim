#!/usr/bin/python

# Script for building plots of output quantities

import argparse
import numpy as np
import matplotlib.pyplot as plt

import pvparse as pvp
import stat_equil as se

parser = argparse.ArgumentParser()

parser.add_argument("file", help="file batch name storing information to plot")
parser.add_argument("-d", "--domain", help="string giving the domain type",
                    required=True)
parser.add_argument("-n", "--nentries", required=True, type=int,
                    help="number of entries in each file")
parser.add_argument("-o", "--out", help="root name for output files")
parser.add_argument("-r", "--nruns", required=True, type=int,
                    help="number of files in the batch")
parser.add_argument("-s", "--scaling", action="store_true",
                    help="generates plots for n, n/2, n/4,... with n=nentries")

args = parser.parse_args()

if (args.out == None):
    root = args.file
else:
    root = args.out
dim = 2
domain = args.domain
nruns = args.nruns
nentries = args.nentries

outfile = args.file + "_%0" + str(len(str(nruns))) + "d" + ".out"

(status, outdict)= pvp.readOutputFileBatch(nruns, nentries, dim, outfile,
                                           compute_mean=False)
# time = np.mean(outdict[pvp.TIME], 0)
time = outdict[pvp.TIME][0]
center = outdict[pvp.CENTER]
moment = outdict[pvp.MOMENT]
hamil = outdict[pvp.HAMIL]
if domain == se.Domain.TWO_LAYER:
    first1 = np.mean(outdict[pvp.CENTER1], 0)
    first2 = np.mean(outdict[pvp.CENTER2], 0)
    second1 = np.mean(outdict[pvp.SECOND1], 0)
    second2 = np.mean(outdict[pvp.SECOND2], 0)
else:
    second = np.mean(outdict[pvp.SECOND], 0)

hamerr = np.zeros((nruns,))
for i in range(nruns):
    hamerr[i] = np.abs(np.max(hamil[i]) - np.min(hamil[i]))
hamil = np.mean(hamil, 0)

# Plot everything
import matplotlib.pyplot as plt
fig, axes = plt.subplots(nrows=dim + 1, ncols=1, figsize=(10,4*(dim + 1)))
for k in range(dim):
    axes[k].plot(time, center[:,:,k].transpose())

axes[dim].plot(center[:,:,0], center[:,:,1], 'bo')

file_str = root + "_center.png"
plt.savefig(file_str)
plt.close()

plt.figure()
plt.plot(time, moment.transpose())
file_str = root + "_moment.png"
plt.savefig(file_str)
plt.close()

plt.figure()
plt.plot(time, hamil)
file_str = root + "_hamil.png"
plt.savefig(file_str)
plt.close()

plt.figure()
plt.plot(range(nruns), hamerr)
file_str = root + "_hamerr.png"
plt.savefig(file_str)
plt.close()

if (domain == se.Domain.TWO_LAYER):
    fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,12))
    axes[0].plot(time, second1[:,0], label="Layer 1")
    axes[0].plot(time, second2[:,0], label="Layer 2")
    axes[0].legend()
    axes[1].plot(time, second1[:,1], label="Layer 1")
    axes[1].plot(time, second2[:,1], label="Layer 2")
    axes[1].legend()
    axes[2].plot(time, second1[:,2], label="Layer 1")
    axes[2].plot(time, second2[:,2], label="Layer 2")
    axes[2].legend()
    fstr = root + "_second_moment.png"
    plt.savefig(fstr, dpi=400)
    plt.close()

    plt.figure()
    plt.plot(time, second1[:,0] - second1[:,2], label="Layer 1")
    plt.plot(time, second2[:,0] - second2[:,2], label="Layer 2")
    plt.plot(time, np.zeros(time.size), 'r--', label="Eq Value")
    plt.legend()
    plt.title("Orientation")
    fstr = root + "_orient.png"
    plt.savefig(fstr, dpi=400)
    plt.close()

    plt.figure()
    plt.plot(time, 2*second1[:,1], label="Layer 1")
    plt.plot(time, 2*second2[:,1], label="Layer 2")
    plt.plot(time, np.zeros(time.size), 'r--', label="Eq Value")
    plt.legend()
    plt.title("Ellipticity")
    fstr = root + "_ellip.png"
    plt.savefig(fstr, dpi=400)
    plt.close()

    file_str = root + "_first_moment.png"
    fig, axes = plt.subplots(nrows=dim, ncols=1, figsize=(10,4*dim))
    for k in range(dim):
        axes[k].plot(time, first1[:,k], label="Layer 1")
        axes[k].plot(time, first2[:,k], label="Layer 2")
        axes[k].legend()
    plt.savefig(file_str, dpi=400)
    plt.close()

else:
    if (args.scaling):
        allsecond = outdict[pvp.SECOND]
        kfstr = "%d Runs"

        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,12))
        n = nruns*2
        while n % 2 == 0:
            n = n / 2
            lstr = kfstr % n
            second = np.mean(allsecond[:n,:,:], 0)
            axes[0].plot(time, second[:,0], label=lstr)
            axes[1].plot(time, second[:,1], label=lstr)
            axes[2].plot(time, second[:,2], label=lstr)

        fstr = root + "_second_moment.png"
        plt.legend()
        plt.savefig(fstr, dpi=400)
        plt.close()

        plt.figure()
        n = nruns*2
        while n % 2 == 0:
            n = n / 2
            lstr = kfstr % n
            second = np.mean(allsecond[:n,:,:], 0)
            plt.plot(time, second[:,0] - second[:,2], label=lstr)

        plt.plot(time, np.zeros(time.size), 'r--', label="Eq Value")
        plt.title("Orientation")
        plt.legend()
        fstr = root + "_orient.png"
        plt.savefig(fstr, dpi=400)
        plt.xlim((200, 500))
        plt.ylim((-0.05, 0.05))
        fstr = root + "_orient_zoom.png"
        plt.savefig(fstr, dpi=400)
        plt.close()

        plt.figure()
        n = nruns*2
        while n % 2 == 0:
            n = n / 2
            lstr = kfstr % n
            second = np.mean(allsecond[:n,:,:], 0)
            plt.plot(time, 2*second[:,1], label=lstr)
        plt.plot(time, np.zeros(time.size), 'r--')
        plt.title("Ellipticity")
        plt.legend()
        fstr = root + "_ellip.png"
        plt.savefig(fstr, dpi=400)
        plt.xlim((200, 500))
        plt.ylim((-0.05, 0.05))
        fstr = root + "_ellip_zoom.png"
        plt.savefig(fstr, dpi=400)
        plt.close()
    else:
        fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,12))
        axes[0].plot(time, second[:,0])
        axes[1].plot(time, second[:,1])
        axes[2].plot(time, second[:,2])
        fstr = root + "_second_moment.png"
        plt.savefig(fstr, dpi=400)
        plt.close()

        plt.figure()
        plt.plot(time, second[:,0] - second[:,2])
        plt.plot(time, np.zeros(time.size), 'r--', label="Eq Value")
        plt.title("x^2 - y^2")
        fstr = root + "_orient.png"
        plt.savefig(fstr, dpi=400)
        plt.close()

        plt.figure()
        plt.plot(time, 2*second[:,1])
        plt.plot(time, np.zeros(time.size), 'r--')
        plt.title("2xy")
        fstr = root + "_ellip.png"
        plt.savefig(fstr, dpi=400)
        plt.close()
