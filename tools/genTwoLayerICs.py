#!/usr/bin/python

# Script for generating initialization files for a point vortex simulation

import argparse
import numpy as np
import sys

import common_observables as co
import optimal_closure as oc
import poisson_solver as ps
import pvparse as pvp
import stat_equil as se

import scipy.integrate
import scipy.interpolate
import scipy.optimize

import matplotlib.pyplot as plt

# BOX = "box"
CIR = "circle"
ENERGY = "energy"
NORM = "norm"
TWO_PATCH = "two-patch"
# UNIF = "unif"

# Functions for building the cdfs
def zeta(y, zeta_interp, x0):
    return zeta_interp(x0, y)
def conditionalZeta(y, zeta, zetax0, x0):
    return zeta(x0, y) / zetax0
def findx(x, zx, p):
    return zx(x) - p
def findy(y, cond_cdf_y, x0, p):
    return cond_cdf_y.ev(x0, y) - p

# Create command line argument parser
parser = argparse.ArgumentParser()

# Define accepted arguments
parser.add_argument("file", help="name of the file to write initial data to")
parser.add_argument("-d", "--dist", default=None, nargs='+', required=True,
                    help="distribution from which to draw initial positions")
parser.add_argument("-n","--number", help="number of point vortices",
                    type=int, required=True)
parser.add_argument("-r", "--nfiles", type=int, default=1,
                    help="number of initial condition files to generate")
parser.add_argument("-s", "--seed", type=int,
                    help="seed for random number generation")
parser.add_argument("-v","--vorticity",
                    help="voriticity of each point vortex", type=float)
parser.add_argument("-l", "--defrad", required=True, type=float,
                    help="deformation radius for the two layer equations")
parser.add_argument("--start", default=0, type=int,
                    help="file number to start from -- use this option for adding to a previous set of runs")
parser.add_argument("-a", "--append", help="append the output to the file",
                    action="store_true")
parser.add_argument("-p", "--plot", help="Plot the initial condition",
                    action="store_true")

# Parse the command line
args = parser.parse_args()

# Build the output vorticity
N = args.number
Nlayer = int(N/2)
if (args.vorticity == None):
    gamma = 2.0/float(N) # Total vorticy 2.0
else:
    gamma = args.vorticity
vort = gamma*np.ones((N,))

if (args.seed != None):
    seed = args.seed
else:
    seed = N

layer = np.ones(N, dtype=np.int)
layer[:Nlayer] *= -1
domain = se.Domain(se.Domain.TWO_LAYER, args.defrad)

# Setup initial distribution stuff
if (args.dist == None):
    distribution = NORM
    cov = np.array([1.0, 0.0, 0.0, 1.0]).reshape((2,2))
    mean = np.zeros((2,))
elif (args.dist[0] == NORM):
    distribution = NORM
    if (len(args.dist) == 1):
        mean = np.zeros((2,))
        cov = np.array([1.0, 0.0, 0.0, 1.0]).reshape((2,2))
    elif (len(args.dist) == 4):
        mean = np.zeros((2,))
        s1 = float(args.dist[1])
        cv = float(args.dist[2])
        s2 = float(args.dist[3])
        cov = np.array([s1, cv, cv, s2]).reshape((2,2))
    elif (len(args.dist) == 6):
        m1 = float(args.dist[1])
        m2 = float(args.dist[2])
        mean = np.array([m1, m2])
        s1 = float(args.dist[3])
        cv = float(args.dist[4])
        s2 = float(args.dist[5])
        cov = np.array([s1, cv, cv, s2]).reshape((2,2))
    else:
        print "--dist norm requires var1 cov var2 or mean1 mean2 var1 cov var2"
        sys.exit(-1)
elif (args.dist[0] == TWO_PATCH):
    distribution = TWO_PATCH
    if (len(args.dist) == 1):
        offset = 0.1
    elif (len(args.dist) == 2):
        offset = float(args.dist[1])
    else:
        print "--dist two-patch [<x-center>]"
        sys.exit(-1)

    # radius = np.sqrt(4.0 - 2.0*offset**2)

    mean1 = np.array([offset, 0.0])
    mean2 = -mean1
    sigsq = (2.0 - offset*offset) / 2.0
    if (sigsq <= 0.0):
        print "Unable to do offset of",offset,"and hold angular impulse at 4.0"
        sys.exit(-1)
    cov1 = np.array([sigsq, 0.0, 0.0, sigsq]).reshape((2,2))
    cov2 = cov1
    grid = ps.Grid(grid_size=100)
    zeta = np.zeros((2,) + grid.xmat.shape)
    zeta[0] = oc.normal_pdf(grid, mean1, cov1)
    zeta[1] = oc.normal_pdf(grid, mean2, cov2)
elif (args.dist[0] == ENERGY):
    distribution = ENERGY

    #### Generate distribution from which to draw samples ####
    print "Computing equilibrium distribution..."

    # Compute the distribution of the initial data
    energy = float(args.dist[1])
    orient = float(args.dist[2])
    ellip = float(args.dist[3])

    grid = ps.Grid(grid_size=64)

    # Build initial vorticity and streamfunction
    cov = np.array([1.0, 0.6, 0.6, 1.0]).reshape((2,2))
    Lsq = cov[0,0] + cov[1,1]
    denom = np.linalg.det(cov)
    beta = 0.0
    alpha = -0.25 * Lsq / (denom)
    mu = -np.log(2.0 * np.pi * np.sqrt(denom))
    lamb1 = 0.25 * (cov[0,0] - cov[1,1]) / denom
    lamb2 = 0.5 * cov[0,1] / denom
    zeta0 = np.exp(alpha * grid.rsq + lamb1 * (grid.xmat**2 - grid.ymat**2)
                   + lamb2 * (2.0 * grid.xmat * grid.ymat) + mu)

    psi0 = ps.solvePoisson(zeta0, grid)
    b0 = -30.0
    a0 = 0.0
    mu1_0 = 0.0
    mu2_0 = 0.0
    z0 = np.zeros((2,) + zeta0.shape)
    z0[0] = z0[1] = zeta0
    p0 = np.zeros((2,) + psi0.shape)
    p0[0] = p0[1] = psi0
    fargs = [grid, p0, z0, a0, b0, mu1_0, mu2_0, energy, 4.0, 1.0, 1.0]

    steq = se.EqDist.buildEqDist(domain, *fargs)
    steq = se.computeStatisticalEquilibrium(steq)
    if steq.ateq == False:
        print "****ERROR: Failed to find statistical equilibrium. Aborting."
        sys.exit(-1)
    else:
        print "Equilibrium Conditions:"
        print "beta =", steq.beta, "alpha =", steq.alpha

    steq.addObservable(0.1, co.A1, orient, name="A1", layer=0)
    steq.addObservable(0.1, co.A2, ellip, name="A2", layer=0)
    steq = se.computeStatisticalEquilibrium(steq)

    if steq.ateq == False:
        print "****ERROR: Failed to find statistical equilibrium. Aborting."
        sys.exit(-1)
    else:
        print "Constrained Conditions:"
        print "beta =", steq.beta, "alpha =", steq.alpha
        print "Done."

    #### Setup densities for generating samples ####
    print "Building necessary marginal and distribution functions..."
    # Density interpolation
    xi = grid.xvect()
    yi = grid.yvect()
    cdf_x = []
    ccdf_y = []
    zeta = steq.zeta

    for k in range(2):
        zeta_interp = scipy.interpolate.RectBivariateSpline(xi, yi,
                                                            steq.zeta[k])

        # Build the cdfs (actually interpolants)
        zeta_x = np.zeros(xi.shape)
        cdf_xa = np.zeros(xi.shape)
        for i in range(xi.size):
            zeta_x[i], err = scipy.integrate.quad(zeta, -np.inf, np.inf,
                                                  args=(zeta_interp, xi[i]),
                                                  limit=100)
        zeta_x = scipy.interpolate.InterpolatedUnivariateSpline(xi, zeta_x,
                                                                ext=1)
        for i in range(xi.size):
            cdf_xa[i], err = scipy.integrate.quad(zeta_x, -np.inf, xi[i],
                                                 limit=100)
        cdf_x.append(
            scipy.interpolate.InterpolatedUnivariateSpline(xi, cdf_xa, ext=1))
        ccdf_ya = np.zeros(grid.xmat.shape)
        for i in range(xi.size):
            for j in range(yi.size):
                val, err = scipy.integrate.quad(conditionalZeta, -np.inf, yi[j],
                                                args=(zeta_interp,
                                                      zeta_x(xi[i]),xi[i]),
                                                limit=200, epsabs=1e-8,
                                                epsrel=1e-8)
                ccdf_ya[i,j] = val
                # if err > 1e-6:
                #     print "Excess error: x =", xi[i], ", y =", yi[j], "error =", err
        ccdf_y.append(scipy.interpolate.RectBivariateSpline(xi, yi, ccdf_ya))
    print "Done."

# Build each file
print "Generating positions..."
root = args.file + "_%0" + str(len(str(args.nfiles))) + "d"

info = "# npvs: " + str(N) + "\n"
info += "# seed: %d\n"
info += "# dist: "
for item in args.dist:
    info += str(item) + " "
info += "\n"
info += "# def radius: " + str(domain.defrad) + "\n"

for k in range(args.start, args.nfiles):
    # Progress bar
    pvp.reportProgress(k - args.start, args.nfiles - args.start, 80)
    # Set seed
    np.random.seed(seed + k)
    # Generate Positions
    if (distribution == NORM):
        pos = np.random.multivariate_normal(mean, cov, N)
    elif (distribution == TWO_PATCH):
        pos = np.zeros((N,2))
        indx1 = np.nonzero(layer == -1)
        indx2 = np.nonzero(layer == 1)
        pos[indx1,:] = np.random.multivariate_normal(mean1, cov1, Nlayer)
        pos[indx2,:] = np.random.multivariate_normal(mean2, cov2, Nlayer)

        # pos = np.zeros((N, 2))
        # rad = np.sqrt(np.random.uniform(low=0.0, high=radius, size=N))
        # ang = np.random.uniform(low=-np.pi, high=np.pi, size=N)
        # pos[:,0] = rad * np.cos(ang)
        # pos[:,1] = rad * np.sin(ang)

    elif (distribution == ENERGY):
        pos = np.zeros((args.number, 2))
        for l in range(2):
            px = np.random.uniform(size=Nlayer)
            py = np.random.uniform(size=Nlayer)
            xcoord = scipy.optimize.fsolve(findx, np.zeros(px.shape),
                                           args=(cdf_x[l], px))
            ycoord = scipy.optimize.fsolve(findy, np.zeros(py.shape),
                                           args=(ccdf_y[l], xcoord, py))
            pos[l*Nlayer:(l+1)*Nlayer,0] = xcoord
            pos[l*Nlayer:(l+1)*Nlayer,1] = ycoord
    else:
        print distribution, "not implemented"
        sys.exit(-1)

    # Force center of vorticity
    if (distribution == TWO_PATCH):
        indx1 = np.nonzero(layer == -1)
        indx2 = np.nonzero(layer ==1)
        # Force center of vorticity to be the desired value in each layer --
        # this forces the origin to be the center of vorticity
        for j in range(2):
            pos[indx1, j] += mean1[j] - np.mean(pos[indx1, j])
            pos[indx2, j] += mean2[j] - np.mean(pos[indx2, j])
    else:
        # Force center of vorticity to be the origin
        for j in range(2):
            pos[:,j] -= np.mean(pos[:,j])

    # Force angular momentum to be 4.0
    lsq = np.sum(vort * np.sum(pos*pos,1))
    ell = np.sqrt(lsq)
    pos *= 2.0/ell
    lsq = np.sum(vort * np.sum(pos*pos,1))

    # Write to file
    fname = root % k
    pvp.writeInitFile(fname, args.append, info % (seed + k), vort, pos, layer)

    # Generate plot if requested
    if (args.plot):
        plt.figure(figsize=(8,8))
        plt.autumn()
        plt.scatter(pos[:Nlayer, 0], pos[:Nlayer, 1], c='b', zorder=2,
                    label="Layer 1")
        plt.scatter(pos[Nlayer:, 0], pos[Nlayer:, 1], c='r', zorder=1,
                    label="Layer 2")
        plt.axes().set_aspect('equal')
        plt.xlim(-4.0, 4.0)
        plt.ylim(-4.0, 4.0)
        plt.legend()
        fstr = fname + ".png"
        plt.savefig(fstr)
        plt.close()

        plt.figure(figsize=(8,8))
        plt.autumn()
        plt.scatter(pos[:Nlayer, 0], pos[:Nlayer, 1], c='b', zorder=1)
        if (distribution == ENERGY or distribution == TWO_PATCH):
            ct = plt.contour(grid.xmat, grid.ymat, zeta[0], linewidths=2.0,
                             zorder=2)
            plt.colorbar(ct, shrink=0.8)
        plt.axes().set_aspect('equal')
        plt.xlim(-4.0, 4.0)
        plt.ylim(-4.0, 4.0)
        fstr = fname + "_L1.png"
        plt.savefig(fstr)
        plt.close()

        plt.figure(figsize=(8,8))
        plt.autumn()
        plt.scatter(pos[Nlayer:, 0], pos[Nlayer:, 1], c='b', zorder=1)
        if (distribution == ENERGY or distribution == TWO_PATCH):
            ct = plt.contour(grid.xmat, grid.ymat, zeta[1], linewidths=2.0,
                             zorder=2)
            plt.colorbar(ct, shrink=0.8)
        plt.axes().set_aspect('equal')
        plt.xlim(-4.0, 4.0)
        plt.ylim(-4.0, 4.0)
        fstr = fname + "_L2.png"
        plt.savefig(fstr)
        plt.close()

pvp.reportProgress(args.nfiles, args.nfiles, 80)
print "Done."
