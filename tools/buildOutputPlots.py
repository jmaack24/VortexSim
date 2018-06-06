#!/usr/bin/python

# Script for building plots of output quantities

import argparse
import numpy as np
import matplotlib.pyplot as plt

import pvparse as pvp

parser = argparse.ArgumentParser()

parser.add_argument("file", help="file storing information to plot")
parser.add_argument("-o", "--out", help="root name for output files")
parser.add_argument("-d", "--dim", help="dimension of quantities", type=int,
                    required=True)
parser.add_argument("-c", "--compare", action="store_true",
                    help="plot covariances on axis with autocorrelation")
parser.add_argument("--start", type=float, help="specify simulation time to start plots at (inclusive)")
parser.add_argument("--stop", type=float,
                    help="specify simulation time to stop plots at (inclusive)")

args = parser.parse_args()

if (args.out == None):
    root, sep, last = args.file.rpartition(".")
else:
    root = args.out
dim = args.dim

status, outdict = pvp.readOutputFile(args.file, dim)

time = outdict[pvp.TIME]
center = outdict[pvp.CENTER]
moment = outdict[pvp.MOMENT]
hamil = outdict[pvp.HAMIL]
second = outdict[pvp.SECOND]
third = outdict[pvp.THIRD]
fourth = outdict[pvp.FOURTH]

if (args.start == None):
    t0 = time[0]
else:
    t0 = args.start
if (args.stop == None):
    tf = time[-1]
else:
    tf = args.stop

# Plot everything
#### Center of Vorticity ####
fig, axes = plt.subplots(nrows=dim + 1, ncols=1, figsize=(10,4*(dim + 1)))
for k in range(dim):
    axes[k].plot(time, center[:,k])
    axes[k].set_xlim([t0, tf])

# TODO Find satisfactory method of only plotting center of vorticity points in
# the given time interval

for k in range(time.size):
    if (time[k] >= t0 and time[k] <= tf):
        axes[dim].plot(center[k,0], center[k,1], 'bo')

# axes[dim].plot(center[:,0], center[:,1], 'bo')

file_str = root + "_center.png"
plt.savefig(file_str)
plt.close()

#### Angular Impluse/Moment of Inertia ####
plt.figure()
plt.plot(time, moment)
plt.xlim([t0, tf])
file_str = root + "_moment.png"
plt.savefig(file_str)
plt.close()

#### Hamiltonian ####
plt.figure()
plt.plot(time, hamil)
plt.xlim([t0, tf])
file_str = root + "_hamil.png"
plt.savefig(file_str)
plt.close()

#### 2nd Moments ####
if (dim == 2):
    fig, axes = plt.subplots(nrows=3, ncols=1)

    axes[0].plot(time, second[:,0])
    axes[0].set_xlim([t0,tf])

    axes[1].plot(time, second[:,1])
    axes[1].set_xlim([t0,tf])

    axes[2].plot(time, second[:,2])
    axes[2].set_xlim([t0,tf])

    fstr = root + "_second_moment.png"
    plt.savefig(fstr)
    plt.close()

#### 3rd Moments ####
if (dim == 2):
    fig, axes = plt.subplots(nrows=4, ncols=1)

    axes[0].plot(time, third[:,0])
    axes[0].set_xlim([t0,tf])

    axes[1].plot(time, third[:,1])
    axes[1].set_xlim([t0,tf])

    axes[2].plot(time, third[:,2])
    axes[2].set_xlim([t0,tf])

    axes[3].plot(time, third[:,3])
    axes[3].set_xlim([t0,tf])

    fstr = root + "_third_moment.png"
    plt.savefig(fstr)
    plt.close()

#### 4th Moments ####
if (dim == 2):
    fig, axes = plt.subplots(nrows=5, ncols=1)

    axes[0].plot(time, fourth[:,0])
    axes[0].set_xlim([t0,tf])

    axes[1].plot(time, fourth[:,1])
    axes[1].set_xlim([t0,tf])

    axes[2].plot(time, fourth[:,2])
    axes[2].set_xlim([t0,tf])

    axes[3].plot(time, fourth[:,3])
    axes[3].set_xlim([t0,tf])

    axes[4].plot(time, fourth[:,4])
    axes[4].set_xlim([t0,tf])

    fstr = root + "_fourth_moment.png"
    plt.savefig(fstr)
    plt.close()

#### Observables ####
if (dim == 2):
    fig, axes = plt.subplots(nrows=2, ncols=1)
    # x^2 - y^2
    axes[0].plot(time, second[:,0] - second[:,2])
    axes[0].title.set_text("x^2 - y^2")
    # 2xy
    axes[1].plot(time, 2*second[:,1])
    axes[1].title.set_text("2xy")

    fstr = root + "_asym.png"
    plt.savefig(fstr)
    plt.close()
