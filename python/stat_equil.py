
import numpy as np
import scipy.optimize

import poisson_solver as ps

#### Useful classes ####

class Domain:
    """Domain and interaction type for point vortex simulation"""

    # Useful quantities
    EBetaZeroPlane = -3.21922184678e-2
    # EBetaZeroScreened = 5.33504013839e-2
    # EBetaZeroTwoLayer = 2.0*EBetaZeroPlane

    # Type names
    PLANE = "plane"
    SCREENED = "screened"
    TWO_LAYER = "two-layer"
    types = [PLANE, SCREENED, TWO_LAYER]

    # Functions
    def __init__(self, name, defrad=1.0):
        self.type = name
        self.dim = 2
        if not name in Domain.types:
            print "Unrecognized domain type::", name
            raise ValueError
        # Type specific information
        if self.type == Domain.SCREENED or self.type == Domain.TWO_LAYER:
            self.defrad = defrad

class ConjugatePair:
    """Coefficient and observable pair for a Gibbs distribution"""
    def __init__(self, coef, obs_func, value, name, eval_func, layer=0):
        # Function to evaluate the observable on a grid object
        self.ofunc = obs_func
        self.efunc = eval_func
        # Value of the observable at equilibrium
        self.value = value
        # Guess for the initial coefficient value -- when at equilibrium,
        # this will contain the correct coefficient value
        self.coef = coef
        # Save the name if given
        self.name = name
        # Layer that this pair applies to (if applicable)
        self.layer = layer
    def __str__(self):
        return self.name + " " + str(self.coef) + " " + str(self.value)
    def obs(self, grid, psi, zeta):
        return self.ofunc(grid, psi, zeta)
    def eval(self, grid, psi, zeta):
        return self.efunc(grid, psi, zeta, self)

class EqDist:
    """Statistical equilibrium distribution for point vortices"""
    @staticmethod
    def buildEqDist(dom, *args):
        retval = None
        if (dom.type == Domain.PLANE or dom.type == Domain.SCREENED):
            retval = EqDist(dom, *args)
        elif(dom.type == Domain.TWO_LAYER):
            retval = TwoLayerEqDist(dom, *args)
        else:
            print "Unimplemented domain type:", dom.type
        return retval

    @staticmethod
    def defaultEval(grid, psi, zeta, conj_pair):
        return np.sum(grid.wmat * conj_pair.obs(grid, psi, zeta) * zeta)

    #### Object Level Functions ####
    def __init__(self, dom, grid, psi, zeta, alpha, beta, mu, Ener,
                 Lsq=2.0, Circ=1.0):
        self.pairs = []
        self.npairs = 0
        self.addObservable(-beta, energy, Ener, name="energy",
                           eval_func=computeHamiltonian)
        self.addObservable(-alpha, angularimpulse, Lsq, name="angular_imp")
        self.addObservable(-mu, circulation, Circ, name="circulation")

        self.grid = grid
        self.domain = dom
        self.psi = psi
        self.zeta = zeta

        self.H = computeHamiltonian(self.grid, self.psi, self.zeta)

        self.ateq = False

        # Shortcuts for conserved quantities
        self.E = Ener
        self.Lsq = Lsq
        self.C = Circ
        self.beta = beta
        self.alpha = alpha
        self.mu = mu

    def addObservable(self, coef, obs_func, value, name=None,
                      eval_func=None, layer=0):
        if eval_func == None:
            eval_func = EqDist.defaultEval
        if name == None:
            name = "Observable" + str(self.npairs)
        cp = ConjugatePair(coef, obs_func, value, name, eval_func)
        self.pairs.append(cp)
        self.npairs += 1
        self.ateq = False
        return self

    def buildZeta(self, coef):
        ex = np.zeros(self.grid.xmat.shape)
        for i in range(self.npairs):
            cp = self.pairs[i]
            ex += coef[i] * cp.obs(self.grid, self.psi, self.zeta)
        return np.exp(ex)

    def buildDistribution(self, coef):
        dist = self.buildZeta(coef)
        dist /= self.C
        return dist

    def coef(self):
        c = np.zeros((self.npairs,))
        for i in range(self.npairs):
            c[i] = self.pairs[i].coef
        return c

    def value(self):
        v = np.zeros((self.npairs,))
        for i in range(self.npairs):
            v[i] = self.pairs[i].value
        return v

    def diagnostic(self):
        print "Equilibriated:", self.ateq
        print "Name   Coef   ExpVal   ActVal"
        err = np.zeros((self.npairs,))
        for i in range(self.npairs):
            fval = self.pairs[i].eval(self.grid, self.psi, self.zeta)
            info = "pairs[" + str(i) + "]:  "
            info += str(self.pairs[i]) + " " + str(fval)
            print info
            val = self.pairs[i].value
            err[i] = fval - val
        print "CE (Constraint Error):", err
        print "|CE|_L^2 =", np.linalg.norm(err)
        print "|CE|_L^inf =", np.max(np.abs(err))

    def equilibriated(self):
        self.beta = -self.pairs[0].coef
        self.alpha = -self.pairs[1].coef
        self.mu = -self.pairs[2].coef
        self.dist = self.zeta / self.C
        self.H = self.pairs[0].eval(self.grid, self.psi, self.zeta)
        # self.Lsq = self.pairs[1].eval(self.grid, self.psi, self.zeta)
        # self.C = self.pairs[2].eval(self.grid, self.psi, self.zeta)
        self.ateq = True
        return self

    def update(self, coef, psi, zeta):
        for i in range(self.npairs):
            self.pairs[i].coef = coef[i]
        self.psi = psi
        self.zeta = zeta
        self.H = self.pairs[0].eval(self.grid, self.psi, self.zeta)
        return self

    def constraint(self, c0):
        f = np.zeros((self.npairs,))
        zeta = self.buildZeta(c0)
        f[0] = (2.0*self.pairs[0].eval(self.grid, self.psi, zeta)
                - self.H - self.E)
        for i in range(1, self.npairs):
            f[i] = (self.pairs[i].eval(self.grid, self.psi, zeta)
                    - self.pairs[i].value)
        return f

    def constraint_jacob(self, c0):
        J = np.zeros((self.npairs, self.npairs))
        zeta = self.buildZeta(c0)
        for i in range(self.npairs):
            f1 = self.pairs[i].obs(self.grid, self.psi, zeta)
            for j in range(self.npairs):
                f2 = self.pairs[j].obs(self.grid, self.psi, zeta)
                J[i,j] = np.sum(self.grid.wmat * f1 * f2 * zeta)
        return J

class TwoLayerEqDist(EqDist):
    @staticmethod
    def defaultEval(grid, psi, zeta, conj_pair):
        return np.sum(grid.wmat * conj_pair.obs(grid, psi, zeta)
                      * zeta[conj_pair.layer])

    def __init__(self, dom, grid, psi, zeta, alpha, beta, mu1, mu2, Ener,
                 Lsq=4.0, Circ1=1.0, Circ2=1.0):
        self.pairs = []
        self.npairs = 0
        self.addObservable(-beta, energy, Ener, name="energy",
                           eval_func=compTwoLayerHamil, layer=-1)
        self.addObservable(-alpha, angularimpulse, Lsq, name="angular_imp",
                           eval_func=compTwoLayerAngImp, layer=-1)
        self.addObservable(-mu1, circulation, Circ1, name="circulation1",
                           layer=0)
        self.addObservable(-mu2, circulation, Circ2, name="circulation2",
                           layer=1)

        self.grid = grid
        self.domain = dom
        self.psi = psi
        self.zeta = zeta

        self.H = compTwoLayerHamil(grid, psi, zeta)

        self.ateq = False

        # Shortcuts for conserved quantities
        self.E = Ener
        self.Lsq = Lsq
        self.C = np.array([Circ1, Circ2])
        self.beta = beta
        self.alpha = alpha
        self.mu1 = mu1
        self.mu2 = mu2

    def addObservable(self, coef, obs_func, value, name=None,
                      eval_func=None, layer=-1):
        if eval_func == None:
            eval_func = TwoLayerEqDist.defaultEval
        if name == None:
            name = "Observable" + str(self.npairs)
        cp = ConjugatePair(coef, obs_func, value, name, eval_func, layer=layer)
        self.pairs.append(cp)
        self.npairs += 1
        self.ateq = False
        return self

    def buildZeta(self, coef):
        zeta = np.zeros((2,) + self.grid.xmat.shape)
        for l in range(2):
            ex = coef[0] * self.C[l] * self.psi[l]
            for i in range(1, self.npairs):
                cp = self.pairs[i]
                if (cp.layer == l):
                    ex += coef[i] * cp.obs(self.grid, self.psi, self.zeta)
                elif (cp.layer == -1):
                    ex += coef[i] * self.C[l] * cp.obs(self.grid, self.psi,
                                                       self.zeta)
            zeta[l] = np.exp(ex)
        return zeta

    def buildDistribution(self, coef):
        dist = self.buildZeta(coef)
        for l in range(2):
            dist[l] /= C[l]
        return dist

    def constraint_jacob(self, c0):
        J = np.zeros((self.npairs, self.npairs))
        # Energy stuff
        J[0,0] = np.sum(self.grid.wmat
                        * (self.C[0]*self.psi[0]**2*self.zeta[0]
                           + self.C[1]*self.psi[1]**2*self.zeta[1]))
        for i in range(1,self.npairs):
            cp = self.pairs[i]
            if cp.layer == -1:
                # Diff energy constraint
                J[0,i] = np.sum(self.grid.wmat
                                * cp.obs(self.grid, self.psi, self.zeta)
                                * (self.C[0]*self.psi[0]*self.zeta[0]
                                   + self.C[1]*self.psi[1]*self.zeta[1]))
                # Diff wrt beta
                J[i,0] = J[0,i]
            else:
                J[0,i] = np.sum(self.grid.wmat
                                * cp.obs(self.grid, self.psi, self.zeta)
                                * self.psi[cp.layer] * self.zeta[cp.layer])
                J[i,0] = np.sum(self.grid.wmat
                                * cp.obs(self.grid, self.psi, self.zeta)
                                * self.C[cp.layer] * self.psi[cp.layer]
                                * self.zeta[cp.layer])
        # Non-energy constraints/coefficients
        for i in range(1,self.npairs):
            cpi = self.pairs[i]
            fi = cpi.obs(self.grid, self.psi, self.zeta)
            for j in range(1,self.npairs):
                cpj = self.pairs[j]
                fj = cpj.obs(self.grid, self.psi, self.zeta)
                if cpi.layer == -1:
                    if cpj.layer == -1:
                        integ = fi * fj * (self.C[0] * self.zeta[0]
                                            + self.C[1] * self.zeta[1])
                    else:
                        integ = fi * fj * self.zeta[cpj.layer]
                else:
                    if cpj.layer == -1:
                        integ = (fi * fj * self.C[cpi.layer]
                                 * self.zeta[cpi.layer])
                    elif cpi.layer == cpj.layer:
                        integ = (fi * fj * self.zeta[cpi.layer])
                    else:
                        integ = 0.0
                J[i,j] = np.sum(self.grid.wmat * integ)
        return J

    def equilibriated(self):
        self.beta = -self.pairs[0].coef
        self.alpha = -self.pairs[1].coef
        self.mu1 = -self.pairs[2].coef
        self.mu2 = -self.pairs[3].coef
        self.dist = np.zeros(self.zeta.shape)
        for l in range(2):
            self.dist[l] = self.zeta[l] / self.C[l]
        self.H = self.pairs[0].eval(self.grid, self.psi, self.zeta)
        self.ateq = True
        return self

# Evaluation functions for the conserved quantities
def energy(grid, psi, zeta):
    return psi
def angularimpulse(grid, psi, zeta):
    return grid.rsq
def circulation(grid, psi, zeta):
    return np.ones(grid.xmat.shape)

#### Common Use Functions ####
def computeHamiltonian(grid, psi, zeta, cp=None):
    H = 0.5 * np.sum(grid.wmat * psi * zeta)
    return H

def computeAngularImpulse(grid, zeta):
    lsq = np.sum(grid.wmat * grid.rsq * zeta)
    return lsq

def computeCirculation(grid, zeta):
    return np.sum(grid.wmat * zeta)

def compTwoLayerHamil(grid, psi, zeta, cp=None):
    H = 0.5 * (np.sum(grid.wmat * np.sum(psi*zeta, axis=0)))
    return H

def compTwoLayerAngImp(grid, psi, zeta, cp=None):
    return np.sum(grid.wmat * grid.rsq * (np.sum(zeta, axis=0)))

def compTwoLayerCirc(grid, zeta):
    C = np.zeros((2,))
    C[0] = np.sum(grid.wmat * zeta[0])
    C[1] = np.sum(grid.wmat * zeta[1])
    return C

### Functions for computing point vortex statistical equilibrium distribution ##

# def dampedNewton(func, jacob, a0, tol=1e-8, niter=1000, args=None):
#     x0 = a0
#     fx0 = func(x0, args)
#     res = np.linalg.norm(fx0, ord=2)
#     count = 0

#     while res > tol:
#         if count > niter:
#             print "***WARNING: Max iterations =", niter, "reached without convergence"
#             break
#         m = 0
#         Jx0 = jacob(x0, args)
#         delta = np.linalg.solve(Jx0, fx0)
#         x1 = x0 - delta
#         fx1 = func(x1, args)
#         err = np.linalg.norm(fx1, ord=2)

#         while err > res:
#             m += 1
#             alpha = 2**(-m)
#             x1 = x0 - alpha*delta
#             fx1 = func(x1, args)
#             err = np.linalg.norm(fx1, ord=2)

#         res = err
#         x0 = x1
#         fx0 = fx1
#         count += 1

#     return x0

def constraint(p0, eq_dist):
    return eq_dist.constraint(p0)

def constraint_jacob(p0, eq_dist):
    return eq_dist.constraint_jacob(p0)

def computeStatisticalEquilibrium(eqdist, rtol=1e-10, mtol=5e-6):
    dom = eqdist.domain
    grid = eqdist.grid

    wmat = grid.wmat
    xmat = grid.xmat
    ymat = grid.ymat
    rsq = grid.rsq

    zeta0 = eqdist.zeta
    psi0 = eqdist.psi
    E0 = eqdist.E
    Lsq = eqdist.Lsq

    # Check that inputs satisfy necessary conditions
    if (dom.type == Domain.TWO_LAYER):
        # H0 = compTwoLayerHamil(grid, psi0, zeta0) # see comment below
        Lsq0 = compTwoLayerAngImp(grid, psi0, zeta0)
        C0 = compTwoLayerCirc(grid, zeta0)
    else:
        # H0 = computeHamiltonian(grid, psi0, zeta0) # see next comment
        Lsq0 = computeAngularImpulse(grid, zeta0)
        C0 = computeCirculation(grid, zeta0)

    # The condition H0 >= E0 is required for the convergence proof. In practice
    # I have found that it does not make a difference and that as long as the
    # initial guess for the coefficients is good enough, the algorithm will
    # converge to the desired equilibrium regardless of the energy of the
    # initial vorticity distribution. As a result I have commented out the
    # check.
    # if (H0 < E0):
    #     print "***WARNING: H(zeta0) >= E0 required"
    #     print "E0 =", E0, ", H(zeta0) =", H0, ", H0 - E0 =", H0 - E0
    #     #raise RuntimeError

    if (np.abs(Lsq0 - Lsq) > 1e-6):
        print "***WARNING: L(zeta0) == Lsq required"
        print "Lsq =", Lsq, " -- L(zeta0) =", Lsq0
        #raise RuntimeError
    if (np.max(np.abs(C0 - eqdist.C)) > 1e-4):
        print "***WARNING: C(zeta0) == 1.0 required"
        print "C(zeta0) =", C0, "    Error =", np.abs(C0 - eqdist.C)
        #raise RuntimeError

    M,N = xmat.shape

    # Setup Hermite differentiation matrix
    if (dom.type == Domain.PLANE):
        Lmat = (-1) * ps.build2dPoissonMat(grid.xvect(), grid.yvect())
        # Get growth solution of Poisson equation
        psiG, zetaG = ps.buildGSoln(grid)
        fp = constraint_jacob
    elif (dom.type == Domain.SCREENED):
        Lmat = ps.build2dScreenedPoissonMat(grid.xvect(), grid.yvect(),
                                            dom.defrad)
        fp = constraint_jacob
    elif (dom.type == Domain.TWO_LAYER):
        Dbarotropic = -ps.build2dPoissonMat(grid.xvect(), grid.yvect())
        Dbaroclinic = ps.build2dScreenedPoissonMat(grid.xvect(), grid.yvect(),
                                                   np.sqrt(0.5) * dom.defrad)
        psi1 = eqdist.psi
        fp = constraint_jacob
    else:
        print "Unimplemented domain/interaction type:", dom.type

    # Setup loop variables
    p1 = eqdist.coef()
    niter = 0
    zeta0_mag = np.max(np.abs(zeta0))
    stop = False #(np.abs((H0 - E0)/E0) <= rtol)
    # Iterate to convergence
    while (not stop):
        # Compute updated parameters

        # print "Updating coefficients..."

        p0 = p1
        p1, info, sts, msg = scipy.optimize.fsolve(constraint, p0,
                                                   fprime=fp,
                                                   full_output=1, maxfev=10000,
                                                   args=(eqdist),
                                                   xtol=1e-12)
        # print "coef =", p1, "-- err =", info["fvec"]

        if sts != 1:
            print "Root finder returned", sts
            print "Message: '", msg
            print "Error:  ", info["fvec"]

        # Compute updated vorticity distribution
        zeta1 = eqdist.buildZeta(p1)

        # print "Computing stream function..."

        # Handle any special processing for the streamfunction calculation
        if (dom.type == Domain.PLANE):
            # Compute updated stream function -- Hermite expansion assumes
            # that the function decays to zero at infinity.  This is untrue here
            # since at infinity psi ~ log r.  Use the linearity of the problem
            # to fix this
            zetaD = zeta1 - zetaG
            psiD = scipy.linalg.solve(Lmat, zetaD.reshape(M*N,)).reshape(M, N)
            psi1 = psiD + psiG
        elif (dom.type == Domain.TWO_LAYER):
            psi1[0], psi1[1] = ps.solveTwoLayerMatrices(zeta1[0], zeta1[1],
                                                        grid, Dbarotropic,
                                                        Dbaroclinic)
        else:
            psi1 = scipy.linalg.solve(Lmat, zeta1.reshape(M*N,)).reshape(M,N)

        # Update the EqDist object
        eqdist.update(p1, psi1, zeta1)

        # Check convergence criteria
        niter += 1
        mag_err = np.max(np.abs(zeta1 - zeta0)) / zeta0_mag
        rel_err = np.abs((eqdist.H - E0)/E0)
        stop = ((rel_err <= rtol and mag_err <= mtol) or
                niter > 100)

        # print "niter =", niter
        # print "H =", eqdist.H, " E0 =", E0
        # print "rel_err =", rel_err
        # print "mag_err =", mag_err

        # Update for iteration
        zeta0 = zeta1
        zeta0_mag = np.max(np.abs(zeta0))

    if niter > 100:
        print "Reached maxiter = 100 without convergence"
        print "H =", eqdist.H, " E0 =", E0
        print "rel_err =", rel_err
        print "mag_err =", mag_err
        eqdist.diagnostic()
    else:
        eqdist.equilibriated()

    return eqdist
