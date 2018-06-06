#!/usr/bin/python

# Script for building plots of observables using optimal closure theory

import argparse
import numpy as np
import matplotlib.pyplot as plt

import common_observables as ob
import optimal_closure as oc
import poisson_solver as ps
import pvparse as pvp
import stat_equil as se

import math
import os
import subprocess

import scipy.optimize

def findExtrema(data):
    maxs = ((data[1:-1] > data[0:-2]) * (data[1:-1] > data[2:]))
    mins = ((data[1:-1] < data[0:-2]) * (data[1:-1] < data[2:]))
    # Prevent any false postive extrema in first 10 seconds of data
    maxs[:10] = False
    mins[:10] = False
    return maxs + mins

def maxExtremaDiff(data1, data2):
    # Finds difference in first peak
    ext1 = findExtrema(data1)
    ext2 = findExtrema(data2)
    # Extrema indexing is shifted one compared to data
    id1 = np.argmax(ext1) + 1
    id2 = np.argmax(ext2) + 1
    diff = np.abs(data1[id1] - data2[id2])
    return diff, id1, id2

GRID_SIZE = 48
WIDTH = 80
pargs = ['b', 'b']

# Set locations for various things
home = os.path.expanduser("~")
SCRIPTS = home + "/Vortex/tools/"
OUTPUT = home + "/Vortex/output/"
SAVE = home + "/Vortex/results/"
LOG = home + "/Vortex/log/"
DATA = home + "/Vortex/data/"
BIN = home + "/Vortex/bin/"
MOVIE = home + "/Vortex/movies/"

# Constant strings
ENER = "energy"
NORM = "norm"
UNIF = "unif"
TWO_PATCH = "two-patch"

parser = argparse.ArgumentParser()

parser.add_argument("root", help="root of output file name to get results from")
parser.add_argument("--single-file", action="store_true",
                   help="flag indicating the file is from a single run -- note that --single-file and -n 1 are not equivalent as the -n 1 option still expects batch processing file name formats")
parser.add_argument("-n", "--nentries", help="number of entries per file",
                    required=True, type=int)
parser.add_argument("-r", "--nfiles", help="number of files in batch",
                    required=True, type=int)
parser.add_argument("-v", "--npvs",
                    help="number of point vortices in simulation",
                    type=int)
parser.add_argument("-d", "--domain", default=None, nargs="+",
                    help="domain and interaction type for PV simulation")
parser.add_argument("--dist", nargs="+", default=None,
                    help="distribution used to initialize vorticity")
parser.add_argument("-m", "--movie", action="store_true",
                    help="create contour movie of trial densities")
parser.add_argument("-o", "--obs", type=int, default=0,
                    help="observable group number")
parser.add_argument("--no-positions", help="just plot trial densities contours",
                    action="store_true")
parser.add_argument("--no-dns", action="store_true",
                    help="just run the closure with no comparison to dns (NOT YET IMPLEMENTED)")

args = parser.parse_args()

nentries = args.nentries
nfiles = args.nfiles
npvs = args.npvs
use_data = False
E0 = None # Energy level of the vorticity distribution -- filled later

if (args.domain == None or args.domain[0] == se.Domain.PLANE):
    domain = se.Domain(se.Domain.PLANE)
elif (args.domain[0] == se.Domain.SCREENED):
    drad = float(args.domain[1])
    domain = se.Domain(se.Domain.SCREENED, defrad=drad)
elif (args.domain[0] == se.Domain.TWO_LAYER):
    drad = float(args.domain[1])
    domain = se.Domain(se.Domain.TWO_LAYER, defrad=drad)
else:
    print args.domain, "is not a recognized type for the domain.  Aborting."
    quit()

if (args.dist == None):
    # Retrieve the data from the first input file -- this is
    # the default behavior
    dfile = DATA + args.root + "/"
    dfile += args.root + "_" + len(str(nfiles)) * "0" + ".pv"
    dist = pvp.getGeneratingDistribution(dfile)
else:
    dist = args.dist

if dist == None:
    quit()

if (dist[0] == NORM):
    distribution = NORM
    if (len(dist) == 1):
        use_data = True
    elif (len(dist) == 4):
        mean = np.zeros((2,))
        s1 = float(dist[1])
        cv = float(dist[2])
        s2 = float(dist[3])
        cov = np.array([s1, cv, cv, s2]).reshape((2,2))
    elif (len(dist) == 6):
        m1 = float(dist[1])
        m2 = float(dist[2])
        mean = np.array([m1, m2])
        s1 = float(dist[3])
        cv = float(dist[4])
        s2 = float(dist[5])
        cov = np.array([s1, cv, cv, s2]).reshape((2,2))
    else:
        print "--dist norm requires var1 cov var2 or mean1 mean2 var1 cov var2"
        quit()
elif (dist[0] == ENER):
    distribution = ENER
    mean = np.zeros((2,))
    cov = np.array([1.0, 0.6, 0.6, 1.0]).reshape((2,2))
    if (len(dist) == 4):
        E0 = float(dist[1])
        use_data = True
    else:
        print "Unexpected formatting for distribution line:\n\t", dist
        quit()
elif (dist[0] == TWO_PATCH):
    distribution = TWO_PATCH
    if (len(dist) == 2):
        offset = float(dist[1])
    else:
        print "--dist two-patch [<x-center>]"
        quit()
    mean1 = np.array([offset, 0.0])
    mean2 = -mean1
    sigsq = (2.0 - offset*offset) / 2.0
    if (sigsq <= 0.0):
        print "Unable to do offset of",offset,"and hold angular impulse at 4.0"
        quit()
    cov1 = np.array([sigsq, 0.0, 0.0, sigsq]).reshape((2,2))
    cov2 = cov1
    zeta_init = oc.two_layer_normal_pdf
    zargs = (1.0, mean1, cov1, mean2, cov2)
else:
    print "Unrecognized distribution --", dist[0]
    quit()

if (args.single_file):
    outfile = args.root
else:
    outfile = OUTPUT + args.root + "/"
    outfile += args.root + "_%0" + str(len(str(nfiles))) + "d.out"

if (args.single_file):
    (status, outdict) = pvp.readOutputFile(outfile, domain.dim)
else:
    (status, outdict) = pvp.readOutputFileBatch(nfiles, nentries, domain.dim,
                                                outfile, compute_mean=True)
time = outdict[pvp.TIME]
center = outdict[pvp.CENTER]
moment = outdict[pvp.MOMENT]
hamil = outdict[pvp.HAMIL]

TIME = time[-1]

if domain.type == se.Domain.TWO_LAYER:
    Lsq = 4.0
else:
    Lsq = 2.0

# Setup initial vorticity distribution stuff
if distribution == NORM or distribution == ENER:
    if domain.type == se.Domain.TWO_LAYER:
        zeta_init = oc.two_layer_normal_pdf
        zargs = (1.0, mean, cov)
    else:
        zeta_init = oc.normal_pdf
        zargs = (mean, cov)
    # if use_data:
    #     # mean = center[0,:]
    #     # if domain.type == se.Domain.TWO_LAYER:
    #     #     cov = np.array([second1[0,0], second1[0,1], second1[0,1],
    #     #                     second1[0,2]]).reshape(2,2)
    #     # else:
    #     #     cov = np.array([second[0,0], second[0,1], second[0,1],
    #     #                     second[0,2]]).reshape(2,2)

    #     # print mean
    #     # print cov

    #     mean = np.array([0.0, 0.0])
    #     cov = np.array([1.0, 0.0, 0.0, 1.0]).reshape(2,2)

# File stuff
if(args.single_file):
    first, sep, last = args.root.rpartition(".")
    root = first
else:
    directory = SAVE + args.root + "/"
    if not os.path.exists(directory):
        os.mkdir(directory)
    root = directory + args.root

#### Setup observables ####
if domain.type == se.Domain.TWO_LAYER:
    second1 = outdict[pvp.SECOND1]
    second2 = outdict[pvp.SECOND2]
    if args.obs == 0:
        obs = np.array([second1[:,0] - second1[:,2], # Layer 1 x^2 - y^2
                        2*second1[:,1] # Layer 1 2xy
        ]).transpose()
        observables = oc.ObservableList(oc.Observable("L1 Orient",
                                                      ob.layer1orient,
                                                      ob.dlay1orient),
                                        oc.Observable("L1 Ellip",
                                                      ob.layer1ellip,
                                                      ob.dlay1ellip))
    elif args.obs == 1:
        obs = np.array([second1[:,0] - second1[:,2], # Layer 1 x^2 - y^2
                        2*second1[:,1], # Layer 1 2xy
                        second2[:,0] - second2[:,2], # Layer 2 x^2 - y^2
                        2*second2[:,1], # Layer 2 2xy
        ]).transpose()
        observables = oc.ObservableList(oc.Observable("L1 Orient",
                                                      ob.layer1orient,
                                                      ob.dlay1orient),
                                        oc.Observable("L1 Ellip",
                                                      ob.layer1ellip,
                                                      ob.dlay1ellip),
                                        oc.Observable("L2 Orient",
                                                      ob.layer2orient,
                                                      ob.dlay2orient),
                                        oc.Observable("L2 Ellip",
                                                      ob.layer2ellip,
                                                      ob.dlay2ellip),
        )
    elif args.obs == 2:
        first1 = outdict[pvp.CENTER1]
        first2 = outdict[pvp.CENTER2]
        obs = np.array([first1[:,0],
                        first1[:,1],
                        first2[:,0],
                        first2[:,1],
        ]).transpose()
        observables = oc.ObservableList(oc.Observable("x1 Center",
                                                      ob.x1m, ob.dx1m),
                                        oc.Observable("y1 Center",
                                                      ob.y1m, ob.dy1m),
                                        oc.Observable("x2 Center",
                                                      ob.x2m, ob.dx2m),
                                        oc.Observable("y2 Center",
                                                      ob.y2m, ob.dy2m),
        )
    elif args.obs == 3:
        first1 = outdict[pvp.CENTER1]
        first2 = outdict[pvp.CENTER2]
        obs = np.array([first1[:,0] - first2[:,0],
                        first1[:,1] - first2[:,1],
                        second1[:,0] - second1[:,2], # Layer 1 x^2 - y^2
                        2*second1[:,1], # Layer 1 2xy
                        second2[:,0] - second2[:,2], # Layer 2 x^2 - y^2
                        2*second2[:,1], # Layer 2 2xy
        ]).transpose()
        observables = oc.ObservableList(oc.Observable("x Moment",
                                                      ob.xjoint, ob.dxjoint),
                                        oc.Observable("y Moment",
                                                      ob.yjoint, ob.dyjoint),
                                        oc.Observable("L1 Orient",
                                                      ob.layer1orient,
                                                      ob.dlay1orient),
                                        oc.Observable("L1 Ellip",
                                                      ob.layer1ellip,
                                                      ob.dlay1ellip),
                                        oc.Observable("L2 Orient",
                                                      ob.layer2orient,
                                                      ob.dlay2orient),
                                        oc.Observable("L2 Ellip",
                                                      ob.layer2ellip,
                                                      ob.dlay2ellip),
        )
    elif args.obs == 4:
        first1 = outdict[pvp.CENTER1]
        first2 = outdict[pvp.CENTER2]
        obs = np.array([first1[:,0] - first2[:,0], # x1 - x2
                        first1[:,1] - first2[:,1], # y1 - y2
        ]).transpose()
        observables = oc.ObservableList(oc.Observable("x Moment",
                                                      ob.xjoint, ob.dxjoint),
                                        oc.Observable("y Moment",
                                                      ob.yjoint, ob.dyjoint),
        )
    else:
        print "Unrecognized observable case number:", args.obs
        quit()
else:
    second = outdict[pvp.SECOND]
    third = outdict[pvp.THIRD]
    fourth = outdict[pvp.FOURTH]
    obs = np.array([second[:,0] - second[:,2], # x^2 - y^2
                    2*second[:,1] # 2xy
    ]).transpose()

    observables = oc.ObservableList(oc.Observable("Orientation", ob.A1, ob.dA1),
                                    oc.Observable("Ellipticity", ob.A2, ob.dA2))

a0 = obs[0]

# Compute the model observable paths
optcl = oc.computeObservablePaths(zeta_init, zargs, domain, a0, observables,
                                  E0=E0, Lsq=Lsq, stop=TIME,
                                  compute_conserved=True,
                                  grid_size=GRID_SIZE)

Eif = optcl.Et
Lif = optcl.Lsqt
Cif = optcl.Circt
aif = optcl.nonstat_obs
stat_eq = optcl.stat_eq
zeta_tilde = optcl.zeta_tilde
psi_tilde = optcl.psi_tilde

# Grid shortcuts
grid = optcl.stat_eq.grid
xmat = grid.xmat
ymat = grid.ymat
rsq = grid.rsq
wmat = grid.wmat

# Print information about energy and angular impulse
if domain.type == se.Domain.TWO_LAYER:
    Eeq = se.compTwoLayerHamil(grid, stat_eq.psi, stat_eq.zeta)
    Lsq_eq = se.compTwoLayerAngImp(grid, stat_eq.psi, stat_eq.zeta)
    Veq = se.compTwoLayerCirc(grid, stat_eq.zeta)
    Ebt = np.zeros((nentries,))
    Ebc = np.zeros((nentries,))
    Ape = np.zeros((nentries,))
    for k in range(nentries):
        Ebt[k] = 0.25 * np.sum(grid.wmat * (zeta_tilde[k,0] + zeta_tilde[k,1])
                               * (psi_tilde[k,0] + psi_tilde[k,1]))
        Ebc[k] = 0.25 * np.sum(grid.wmat * (zeta_tilde[k,0] - zeta_tilde[k,1])
                               * (psi_tilde[k,0] - psi_tilde[k,1]))
        Ape[k] = np.sum(grid.wmat * (psi_tilde[k,0] - psi_tilde[k,1])**2)
    Ape *= 0.5 / domain.defrad**2

    # Compute d/dt APE
    PE = np.zeros((observables.size, observables.size))
    alpha = optcl.X[1,:]
    beta = optcl.X[0,:]
    psi_d = optcl.psi_d
    psi_da = optcl.psi_da
    psi_db = optcl.psi_db
    for k in range(observables.size):
        psi_k = psi_d[k] + alpha[k]*psi_da + beta[k]*psi_db

        # temp = np.sum(wmat * psi_k[0] * stat_eq.psi)
        # print temp
        # temp = np.sum(wmat * psi_k[1] * stat_eq.psi)
        # print temp
        # temp = np.sum(wmat * (psi_k[0] - psi_k[1]) * stat_eq.psi)
        # print temp

        for l in range(observables.size):
            psi_l = psi_d[l] + alpha[l]*psi_da + beta[l]*psi_db
            PE[k,l] = np.sum(wmat * (psi_k[0] - psi_k[1])
                             * (psi_l[0] - psi_l[1]))
    PE /= domain.defrad**2
    Ape_dt = np.zeros((nentries,))
    Cinv = optcl.Cinv
    J = optcl.J
    D = optcl.D
    for k in range(nentries):
        M = oc.computeM(time[k], Cinv, J, D)
        lamb = optcl.lambs[k]
        lambdot = np.dot(Cinv, np.dot((J - M), lamb))
        Ape_dt[k] = np.dot(lamb, np.dot(PE, lambdot))

    # # Check for BCE file and parse the info found
    # bcefile = OUTPUT + "/" + args.root + "/" + args.root + ".bce"
    # if (os.path.exists(bcefile)):
    #     sts, temp, bce = pvp.parseBceFile(bcefile, nentries)

    # Check for APE file and parse the info found
    apefile = OUTPUT + "/" + args.root + "/" + args.root + ".ape"
    if (os.path.exists(apefile)):
        sts, temp, ape = pvp.parseApeFile(apefile, nentries)

else:
    Eeq = se.computeHamiltonian(grid, stat_eq.psi, stat_eq.zeta)
    Lsq_eq = se.computeAngularImpulse(grid, stat_eq.zeta)
    Veq = se.computeCirculation(grid, stat_eq.zeta)

max_coef = np.max(np.abs(optcl.lambs))
(tau, idx) = np.unravel_index(np.argmax(np.abs(optcl.lambs)),
                              dims=optcl.lambs.shape)
tau = time[tau]

mag_lambs = np.linalg.norm(optcl.lambs, axis=1)
max_l2c = np.max(mag_lambs)
idx2 = np.argmax(mag_lambs)
tau2 = time[idx2]

M = oc.computeM(TIME, optcl.Cinv, optcl.J, optcl.D)
total_res = np.dot(optcl.lambs[0], np.dot(M, optcl.lambs[0]))

info = "Equil Energy: " + str(Eeq) + "  Equilibrium AI: " + str(Lsq_eq)
info += "  Equil Vorticity: " + str(Veq) + "\n"
info += "beta = " + str(stat_eq.beta) + "  alpha = " + str(stat_eq.alpha)
if domain.type == se.Domain.TWO_LAYER:
    info += "  mu1 = " + str(stat_eq.mu1) + "  mu2 = " + str(stat_eq.mu2) + "\n"
else:
    info += "  mu = " + str(stat_eq.mu) + "\n"
info += "Initial PV Hamiltonian: " + str(hamil[0]) + "\n"
info += "Final PV Hamiltonian:   " + str(hamil[-1]) + "\n"
info += "Max_{t,i} |lambda_i(t)| = " + str(max_coef)
info += "  t = " + str(tau) + "  i = " + str(idx) + "\n"
info += "Max_{t} |lambda(t)|_l2 = " + str(max_l2c)
info += "  t = " + str(tau2) + "  lambda = " + str(optcl.lambs[idx2]) + "\n"
info += "v(lambda, T) = " + str(total_res) + "  T = " + str(TIME) + "\n"
if domain.type == se.Domain.TWO_LAYER:
    apem = np.max(np.abs(Ape_dt))
    info += "Max_t |d/dt APE| = " + str(apem)
    info += "  t = " + str(time[np.argmax(np.abs(Ape_dt))])
    info += "  t_scale = " + str(np.max(Ape)/apem) + "\n"

# Energy info
emin = np.min(Eif)
emax = np.max(Eif)
info += "Emin = " + str(emin) + "  Emax = " + str(emax)
info += "  Diff = " + str(emax - emin)
info += "  Rel Diff = " + str(np.abs((emax - emin)/Eeq)) + "\n"

# Angular impluse info
lmin = np.min(Lif)
lmax = np.max(Lif)
info += "Lmin = " + str(lmin) + "  Lmax = " + str(lmax)
info += "  Diff = " + str(lmax - Lsq)
info += "  Rel Diff = " + str(np.abs((lmax - Lsq)/Lsq)) + "\n"

# Circulation info
cmin = np.min(Cif)
cmax = np.max(Cif)
info += "Cmin = " + str(cmin) + "  Cmax = " + str(cmax)
info += "  Diff = " + str(cmax - cmin)
info += "  Rel Diff = " + str(np.abs((cmin - cmax))/Veq) + "\n"

#errStat = np.abs(obs - ais)
err = np.abs(obs - aif)
for i in range(a0.size):
    denom = np.max(np.abs(obs[:,i]))
    maxerr = np.max(np.abs(err[:,i]))
    if domain.type == se.Domain.TWO_LAYER:
        obs_var = np.sum(wmat * np.sum(observables.eval(i,grid)**2 *
                                       stat_eq.dist, axis=0))
    else:
        obs_var = np.sum(wmat * observables.eval(i, grid)**2 * stat_eq.dist)
    perr = np.sqrt(2.0 * obs_var * total_res)

    info += str(observables[i])
    info += " max = " + str(denom) + " Max Err = " + str(maxerr)
    info += " Rel Err = " + str(maxerr/denom)
    info += " Pred Err = " + str(perr) + "\n"

print info

info += "\nTable Info\n"
info += "Observable  Peak_Obs  Peak_Pred  Abs_Peak_Err  Rel_Peak_Err  "
info += "Max_Abs_Err  Rel_Err\n"
fmt = "$A_{:1d}$ & {:>.5e} & {:>.5e} & {:>.5e} & {:>.5e} & {:>.5e} & {:>.5e}"

for i in range(a0.size):
    denom = np.max(np.abs(obs[:,i]))
    perr, j, k = maxExtremaDiff(obs[:,i], aif[:,i])
    otruth = np.abs(obs[j,i])
    ptruth = np.abs(aif[k,i])
    relperr = perr/otruth
    maxerr = np.max(np.abs(err[:,i]))
    relerr = maxerr/denom
    info += fmt.format(i+1, otruth, ptruth, perr, relperr, maxerr, relerr)
    info += " \\\\\n"

ifile = open(root + "_info.txt", mode="w")
ifile.write(info)
ifile.flush()
ifile.close()

plt.figure()
plt.plot(time, Lsq*np.ones(moment.shape))
plt.plot(time, moment)
plt.plot(time, Lif)
plt.savefig(root + "_moment.png")
plt.close()

plt.figure()
plt.plot(time, Eeq*np.ones(time.size))
plt.plot(time, hamil)
plt.plot(time, Eif)
plt.savefig(root + "_energy.png")
plt.close()

if domain.type == se.Domain.TWO_LAYER:
    fig, axes = plt.subplots(nrows=2, ncols=1)
    axes[0].plot(time, Ebt, 'b')
    axes[0].plot(time, Eif, 'g')
    axes[0].title.set_text("Barotropic")
    axes[1].plot(time, Ebc, 'm')
    axes[1].title.set_text("Baroclinic")
    #axes[1].plot(time, Eif - Eeq, 'r--')
    #axes[1].plot(time, Ebt_dns, 'b')
    plt.tight_layout()
    plt.savefig(root + "_partition_energy.png")
    plt.close()

    fig, axes = plt.subplots(nrows=2, ncols=1)
    # axes[0].plot(time, Ebc, 'm', label="Total")
    if (os.path.exists(apefile)):
        axes[0].plot(time, ape, 'b', label="EDNS")
        ape_rate = (ape[1:] - ape[:-1]) / (time[1:] - time[:-1])
        axes[1].plot(0.5 * (time[1:] - time[:-1]) + time[1:],
                     ape_rate, 'b', label="EDNS")

    axes[0].plot(time, Ape, 'r', label="Closure")
    axes[0].legend()
    axes[0].title.set_text("Available Potential Energy")
    axes[1].plot(time, Ape_dt, 'r', label="Closure")
    axes[1].legend()
    axes[1].title.set_text("APE Rate")
    plt.tight_layout()
    plt.savefig(root + "_ape.png", dpi=400)
    plt.close()


# r4eq = np.sum(grid.wmat * grid.rsq**2 * stat_eq.zeta)
# plt.figure()
# plt.plot(time, r4eq*np.ones(time.shape), 'r--')
# plt.plot(time, fourth[:,0] + 2.0*fourth[:,2] + fourth[:,4], 'b-')
# # plt.plot(time, aif[:,2] + 2.0*aif[:,4] + aif[:,6], 'g-')
# plt.savefig(root + "_r4.png")
# plt.close()

for i in range(a0.size):
    if domain.type == se.Domain.TWO_LAYER:
        eqval = np.sum(np.tile(grid.wmat, (2,1,1)) *
                       observables.eval(i, grid) * stat_eq.zeta)
    else:
        eqval = np.sum(grid.wmat * observables.eval(i, grid) * stat_eq.zeta)

    plt.figure()
    plt.plot(time, obs[:,i], label="EDNS")
    # plt.plot(time, lstsqs[:,i], "r-")
    plt.plot(time, aif[:,i], "g-", label="Closure")
    plt.plot(time, eqval*np.ones(time.shape), "r--", label="Eq Value")
    plt.legend()
    plt.title(str(observables[i]))
    fstr = root + "_" + str(observables[i]) + ".png"
    plt.savefig(fstr, dpi=400)
    plt.close()

    fig, axes = plt.subplots(nrows=2, ncols=1)
    axes[0].plot(time, obs[:,i], label="EDNS")
    #axes[0].plot(time, ais[:,i], "m-")
    axes[0].plot(time, aif[:,i], "g-", label="Closure")
    axes[0].plot(time, eqval*np.ones(time.shape), "r--", label="Eq Value")
    axes[0].title.set_text(str(observables[i]))
    axes[0].legend()
    #axes[1].plot(time, errStat[:,i], "m-")
    axes[1].plot(time, err[:,i], "g-")
    axes[1].title.set_text("Error")

    fstr = root + "_" + str(observables[i]) + "_err.png"
    plt.savefig(fstr)
    plt.close()

if args.movie:
    print "Creating contour plots..."
    # Create scale for colors
    max_vort = np.max(zeta_tilde)
    dx = 0.1 * max_vort
    levels=np.arange(0.0, max_vort + dx, dx)
    # File numbering scheme that will be appended to root name
    format_str = "_dist_%0" + str(len(str(args.nentries))) + "d.png"

    # Open position data file
    if not args.no_positions:
        data_file = MOVIE + "/" + args.root + "/" + args.root + "_"
        data_file += "0"*len(str(nfiles)) + ".pos"
        pos = open(data_file, "r")
        fne = pvp.findNumEntries(pos)
        if (fne != nentries):
            print "Number of positions in data file", data_file, "has", fne, "entries. Expected", nentries, "Unable to make movie."
            quit()

    for k in range(time.size):
        pvp.reportProgress(k, time.size, WIDTH)

        # File specific stuff
        fstr = format_str % k
        title = "Simulation Time:%8.3f" % time[k]

        if not args.no_positions:
            # Get position plots of point vortices
            (sts, num, dim, tm, ly, vort, position) = pvp.readPositionEntry(pos)
            if (sts != 0):
                print "**ERROR:: readPositionEntry returned", sts
                quit()

        if domain.type == se.Domain.TWO_LAYER:
            fig, axes = plt.subplots(nrows=1, ncols=2,
                                     figsize=(8*2, 8))
            plt.autumn()
            for l in range(2):
                if not args.no_positions:
                    lys = np.unique(ly)
                    indx = np.nonzero(ly == pvp.INV_LAYER_MAP[l + 1])
                    axes[l].scatter(position[indx,0], position[indx,1],
                                    c=pargs[l])
                axes[l].set_title('Layer ' + str(l + 1))
                ct = axes[l].contour(xmat, ymat, zeta_tilde[k,l], levels=levels,
                                     linewidths=2.0, zorder=2)
                axes[l].set_xlim(-4.0, 4.0)
                axes[l].set_ylim(-4.0, 4.0)
                axes[l].set_aspect('equal')

            fig.subplots_adjust(right=0.85)
            cbar_ax = fig.add_axes([0.9, 0.15, 0.025, 0.7])
            fig.colorbar(ct, cax=cbar_ax, ticks=levels)
        else:
            # Plot setup stuff
            plt.figure()
            plt.autumn()

            if not args.no_positions:
                # Plot positions
                plt.scatter(position[:,0], position[:,1], zorder=1, s=10, c='b')

            # Plot trial density contours
            ct = plt.contour(xmat, ymat, zeta_tilde[k], levels=levels,
                             linewidths=2.0, zorder=2)
            plt.colorbar(ct, shrink=0.8, ticks=levels)
            plt.xlim(-4.0, 4.0)
            plt.ylim(-4.0, 4.0)

        plt.suptitle(title)
        plt.savefig(root + fstr, dpi=400)
        plt.close()

    pvp.reportProgress(time.size, time.size, WIDTH)

    print "Building movie..."
    files = root + format_str
    movie = root
    if args.no_positions:
        movie += "_dist"
    movie += ".mp4"
    if os.path.exists(movie):
        os.remove(movie)
    pargs = ["ffmpeg", "-i", files, "-c:v", "libx264", "-pix_fmt",
             "yuv420p", movie]
    subprocess.call(pargs)
    # Remove image files
    for k in range(time.size):
        if k % 50 != 0:
            file_str = root + format_str % k
            os.remove(file_str)

    print "Done!"
