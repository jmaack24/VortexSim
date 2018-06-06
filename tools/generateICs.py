#!/usr/bin/python

# Script for generating initialization files for a point vortex simulation

import argparse
import numpy as np
import sys

import common_observables as co
import poisson_solver as ps
import pvparse as pvp
import stat_equil as se

import scipy.integrate
import scipy.interpolate
import scipy.optimize

import matplotlib.pyplot as plt

BOX = "box"
CIR = "circle"
ENERGY = "energy"
NORM = "norm"
TWO_PATCH = "two-patch"
UNIF = "unif"

B0ENERGY = -3.21922184678e-2

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
parser.add_argument("--domain", default=[se.Domain.PLANE], nargs='+',
                    help="domain of the desired equilibrium")
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
if (args.vorticity == None):
    gamma = 1.0/float(N)
else:
    gamma = args.vorticity
vort = gamma*np.ones((N,))

if (args.seed != None):
    seed = args.seed
else:
    seed = N

if args.domain[0] == se.Domain.SCREENED:
    defrad = float(args.domain[1])
    domain = se.Domain(args.domain[0], defrad)
    layer = None
else:
    layer = None
    domain = se.Domain(args.domain[0])

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
elif (args.dist[0] == UNIF):
    distribution = UNIF
    if (len(args.dist) == 1):
        xlim = 2.0
        ylim = 2.0
    elif (args.dist[1] == CIR):
        distribution += CIR
        radius = float(args.dist[2])
    elif (args.dist[1] == BOX):
        distribution += BOX
        xlim = float(args.dist[2])
        ylim = float(args.dist[3])
    else:
        print "--dist unif [box/circle] [<x-axis spread> <y-axis spread>]/[radius]"
        sys.exit(-1)
elif (args.dist[0] == TWO_PATCH):
    distribution = TWO_PATCH
    if (len(args.dist) == 1):
        mean1 = np.array([-0.5, 0.0])
        mean2 = -mean1
        cov1 = np.array([1.0, 0.0, 0.0, 1.0]).reshape((2,2))
        cov2 = cov1
    elif (len(args.dist) == 6):
        mean1 = np.array([float(args.dist[1]), float(args.dist[2])])
        cov1 = np.array([float(args.dist[3]), float(args.dist[4]),
                         float(args.dist[4]), float(args.dist[5])]
        ).reshape((2,2))
        mean2 = -mean1
        cov2 = cov1
    else:
        print "--dist two-patch [<x-center> <y-center> <x-variance> <xy-covariance> <y-variance>]"
        sys.exit(-1)
elif (args.dist[0] == ENERGY):
    distribution = ENERGY

    #### Generate distribution from which to draw samples ####
    print "Computing equilibrium distribution..."

    # Compute the distribution of the initial data
    energy = float(args.dist[1])
    orient = float(args.dist[2])
    ellip = float(args.dist[3])

    grid = ps.Grid(grid_size=48)

    # Build initial vorticity and streamfunction
    cov = np.array([1.0, 0.0, 0.0, 1.0]).reshape((2,2))
    Lsq = cov[0,0] + cov[1,1]
    denom = np.linalg.det(cov)
    beta = 0.0
    alpha = -0.25 * Lsq / (denom)
    mu = -np.log(2.0 * np.pi * np.sqrt(denom))
    lamb1 = 0.25 * (cov[0,0] - cov[1,1]) / denom
    lamb2 = 0.5 * cov[0,1] / denom
    zeta0 = np.exp(alpha * grid.rsq + lamb1 * (grid.xmat**2 - grid.ymat**2)
                   + lamb2 * (2.0 * grid.xmat * grid.ymat) + mu)

    if domain.type == se.Domain.PLANE:
        psi0 = ps.solvePoisson(zeta0, grid)
        b0 = beta
        a0 = alpha
        mu0 = mu
        fargs = [grid, psi0, zeta0, a0, b0, mu0, energy]
    elif domain.type == se.Domain.SCREENED:
        psi0 = ps.solveScreenedPoisson(zeta0, grid, domain.defrad)
        b0 = beta
        a0 = alpha
        mu0 = mu
        fargs = [grid, psi0, zeta0, a0, b0, mu0, energy]

    steq = se.EqDist.buildEqDist(domain, *fargs)
    steq = se.computeStatisticalEquilibrium(steq)
    steq.addObservable(0.1, co.A1, orient, name="A1")
    steq.addObservable(0.1, co.A2, ellip, name="A2")
    steq = se.computeStatisticalEquilibrium(steq, rtol=1e-10, mtol=1e-5)

    if steq.ateq == False:
        print "****ERROR: Failed to find statistical equilibrium. Aborting."
        sys.exit(-1)
    else:
        print "beta =", steq.beta, "alpha =", steq.alpha
        print "Done."

    #### Setup densities for generating samples ####
    print "Building necessary marginal and distribution functions..."
    # Density interpolation
    xi = grid.xvect()
    yi = grid.yvect()
    zeta_interp = scipy.interpolate.RectBivariateSpline(xi, yi, steq.zeta)

    # Build the cdfs (actually interpolants)
    zeta_x = np.zeros(xi.shape)
    cdf_x = np.zeros(xi.shape)
    print "Building x-marginal..."
    for i in range(xi.size):
        zeta_x[i], err = scipy.integrate.quad(zeta, -np.inf, np.inf,
                                              args=(zeta_interp, xi[i]),
                                              limit=100)
    zeta_x = scipy.interpolate.InterpolatedUnivariateSpline(xi, zeta_x, ext=1)

    print "Building x-cdf..."
    for i in range(xi.size):
        cdf_x[i], err = scipy.integrate.quad(zeta_x, -np.inf, xi[i], limit=100)

    # plt.figure(figsize=(8,8))
    # plt.plot(grid.xvect(), cdf_x)
    # plt.show()
    cdf_x = scipy.interpolate.InterpolatedUnivariateSpline(xi, cdf_x, ext=1)
    ccdf_y = np.zeros(grid.xmat.shape)

    print "Building cdf..."
    for i in range(xi.size):
        for j in range(yi.size):
            #print "i =", i, "j =", j
            val, err = scipy.integrate.quad(conditionalZeta, -np.inf, yi[j],
                                            args=(zeta_interp,
                                                  zeta_x(xi[i]),xi[i]),
                                            limit=500, epsabs=1e-8, epsrel=1e-8)
            if val < 0.0:
                #print "value =", val, "< 0 -- replaced by 0"
                val = 0.0
            if val > 1.0:
                #print "value =", val, "> 1 -- replaced by 1"
                val = 1.0
            ccdf_y[i,j] = val
            # if err > 1e-6:
            #     print "Excess error: x =", xi[i], ", y =", yi[j], "error =", err

    # plt.figure(figsize=(8,8))
    # ct = plt.contour(grid.xmat, grid.ymat, ccdf_y, linewidths=2.0,
    #                  zorder=2)
    # plt.colorbar(ct, shrink=0.8)
    # plt.show()

    ccdf_y = scipy.interpolate.RectBivariateSpline(xi, yi, ccdf_y)
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

for k in range(args.start, args.nfiles):
    # Progress bar
    pvp.reportProgress(k - args.start, args.nfiles - args.start, 80)
    # Set seed
    np.random.seed(seed + k)
    # Generate Positions
    if (distribution == NORM):
        pos = np.random.multivariate_normal(mean, cov, N)
    elif (distribution == UNIF + BOX):
        pos = np.zeros((N, dim))
        pos[:,0] = np.random.uniform(low=-xlim, high=xlim, size=N)
        pos[:,1] = np.random.uniform(low=-ylim, high=ylim, size=N)
    elif (distribution == UNIF + CIR):
        pos = np.zeros((N, dim))
        rad = np.sqrt(np.random.uniform(low=0.0, high=radius, size=N))
        ang = np.random.uniform(low=-np.pi, high=np.pi, size=N)
        pos[:,0] = rad * np.cos(ang)
        pos[:,1] = rad * np.sin(ang)
    elif (distribution == TWO_PATCH):
        pos = np.zeros((N,dim))
        pos[:N/2,:] = np.random.multivariate_normal(mean1, cov1, N/2)
        pos[N/2:,:] = np.random.multivariate_normal(mean2, cov2, N - N/2)
    elif (distribution == ENERGY):
        pos = np.zeros((args.number, 2))
        px = np.random.uniform(size=args.number)
        py = np.random.uniform(size=args.number)
        xcoord = scipy.optimize.fsolve(findx, np.zeros(px.shape),
                                       args=(cdf_x, px))
        ycoord = scipy.optimize.fsolve(findy, np.zeros(py.shape),
                                       args=(ccdf_y, xcoord, py))
        if np.max(np.abs(xcoord)) > 5.0 or np.max(np.abs(ycoord)) > 5.0:
            print np.max(np.abs(xcoord)), np.max(np.abs(ycoord))
        pos[:,0] = xcoord
        pos[:,1] = ycoord
    else:
        print distribution, "not implemented"
        sys.exit(-1)

    # Force center of vorticity to be the origin
    for j in range(2):
        pos[:,j] -= np.mean(pos[:,j])

    # Force angular momentum to be 2.0
    lsq = np.sum(vort * np.sum(pos*pos,1))
    ell = np.sqrt(lsq)
    pos *= np.sqrt(2.0)/ell

    # Write to file
    fname = root % k
    pvp.writeInitFile(fname, args.append, info % (seed + k), vort, pos, layer)

    # Generate plot if requested
    if (args.plot):
        plt.figure(figsize=(8,8))
        plt.autumn()
        plt.scatter(pos[:,0], pos[:,1], c='b', zorder=1)
        if (distribution == ENERGY):
            # Plot actual density contour
            ct = plt.contour(grid.xmat, grid.ymat, steq.zeta, linewidths=2.0,
                             zorder=2)
            plt.colorbar(ct, shrink=0.8)
        plt.xlim(-5.0, 5.0)
        plt.ylim(-5.0, 5.0)
        fstr = fname + ".png"
        plt.savefig(fstr)
        plt.close()

pvp.reportProgress(args.nfiles, args.nfiles, 80)
print "Done."
