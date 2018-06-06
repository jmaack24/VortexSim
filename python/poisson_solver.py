
import numpy as np

#### Helper class ####

class Grid:
    def __init__(self, grid_size=48):
        # Get points and weights for gauss hermite quadrature
        x, wx = np.polynomial.hermite.hermgauss(grid_size)
        y, wy = np.polynomial.hermite.hermgauss(grid_size)

        # Save 1d weights
        self.wx = wx
        self.wy = wy

        # Setup grid
        self.xmat = np.repeat(x, y.size, axis=0).reshape(x.size, y.size)
        self.ymat = np.repeat(y, x.size, axis=0).reshape(y.size,
                                                         x.size).transpose()
        self.rsq = self.xmat*self.xmat + self.ymat*self.ymat

        self.wmat = (np.repeat(wx, y.size, axis=0).reshape(x.size, y.size) *
                     np.repeat(wy, x.size, axis=0).reshape(y.size,
                                                           x.size).transpose()*
                     np.exp(self.rsq))
    def xvect(self):
        return self.xmat[:,0]
    def yvect(self):
        return self.ymat[0,:]

class LegGrid:
    def __init__(self, grid_size=48, box=1.0):
        # Get points and weights for gauss hermite quadrature
        x, wx = np.polynomial.legendre.leggauss(grid_size)
        y, wy = np.polynomial.legendre.leggauss(grid_size)

        wx *= box
        wy *= box
        x *= box
        y *= box

        # Save 1d weights
        self.wx = wx
        self.wy = wy

        # Setup grid
        self.xmat = np.repeat(x, y.size, axis=0).reshape(x.size, y.size)
        self.ymat = np.repeat(y, x.size, axis=0).reshape(y.size,
                                                         x.size).transpose()
        self.rsq = self.xmat*self.xmat + self.ymat*self.ymat

        self.wmat = (np.repeat(wx, y.size, axis=0).reshape(x.size, y.size) *
                     np.repeat(wy, x.size, axis=0).reshape(y.size,
                                                           x.size).transpose())
    def xvect(self):
        return self.xmat[:,0]
    def yvect(self):
        return self.ymat[0,:]

class TwoLayerGrid(Grid):
    def __init__(self, grid_size=48):
        super(TwoLayerGrid, self).__init__(grid_size)
        self.wmat = np.tile(self.wmat, (2,1,1))
        self.xmat = np.tile(self.xmat, (2,1,1))
        self.ymat = np.tile(self.ymat, (2,1,1))
        self.rsq = np.tile(self.rsq, (2,1,1))
    def xvect(self):
        return self.xmat[0,:,0]
    def yvect(self):
        return self.ymat[0,0,:]

#### Functions for solving Poisson's equation ####

def hermdm(xi, nder):
    npts = int(xi.size)
    rettuple = (xi,)

    # Hermite specific stuff
    alpha = np.exp(-xi*xi/2)
    beta = np.ones((nder+1, npts))
    beta[1,:] = -xi
    for l in range(2,nder+1):
        beta[l,:] = -xi*beta[l-1,:] - (l - 1)*beta[l-2,:]

    beta = np.delete(beta, 0, 0)

    # General differentiation matrix stuff
    xmat = np.repeat(xi, npts, axis=0).reshape(npts,npts)
    dx = xmat - xmat.transpose()
    np.fill_diagonal(dx, 1)
    c = alpha * np.prod(dx, 1)
    cmat = np.repeat(c, npts, axis=0).reshape(npts,npts)
    cmat = cmat / cmat.transpose()
    zmat = 1.0/dx
    np.fill_diagonal(zmat, 0)

    D = np.eye(npts)
    ymat = np.ones((npts, npts))
    Xmat = np.zeros((npts - 1, npts))
    for j in range(npts):
        Xmat[:j,j] = zmat[j,:j]
        Xmat[j:,j] = zmat[j,j+1:]

    for l in range(1, nder+1):
        DD = np.repeat(np.diagonal(D), npts, axis=0).reshape(npts,npts)
        D = l * zmat * (cmat*DD - D)
        #np.fill_diagonal(D, 0)
        ymat = np.cumsum(np.vstack(
            ([beta[l-1,:], l*ymat[0:npts-1,:]*Xmat])),axis=0)
        D = D + np.diag(ymat[-1,:])
        rettuple += (D,)
    return rettuple

def build2dPoissonMat(x, y):
    M = x.size
    N = y.size
    xs, D1x, D2x = hermdm(x, 2)
    ys, D1y, D2y = hermdm(y, 2)
    L = np.kron(D2x, np.eye(N)) + np.kron(np.eye(M), D2y)
    return L

def build2dScreenedPoissonMat(x, y, defrad):
    L = build2dPoissonMat(x, y)
    S = -L + 1.0/(defrad**2) * np.eye(x.size * y.size)
    return S

def build2dTwoLayerPoissonMat(x, y, defrad):
    M = x.size
    N = y.size
    f = 1.0 / defrad
    L = np.zeros((2 * M * N, 2 * M * N))
    L1 = L2 = -1.0 * build2dPoissonMat(x, y)
    L[:M*N, :M*N] = L1 + f**2 * np.eye(M * N)
    L[:M*N, M*N:] = -f**2 * np.eye(M * N)
    L[M*N:, M*N:] = L2 + f**2 * np.eye(M * N)
    L[M*N:, :M*N] = -f**2 * np.eye(M * N)
    return L

def buildGSoln(grid, l=5.5):
    # Builds a smooth solution of -Delta psi = zeta which looks like
    #    -1/(2*pi) * log r
    # for r large
    lsq = l*l
    etasq = grid.rsq/(lsq)
    psi = -0.25/np.pi*np.log(grid.rsq)*(1-np.exp(-etasq*etasq))
    zeta = (4.0/(lsq*np.pi)*np.exp(-etasq*etasq)*
            etasq*(1 + (1 - etasq*etasq)*np.log(grid.rsq)))
    return (psi, zeta)

def buildGFlow(grid, l=5.5):
    # Builds psi'(r) * 1/r where psi is from the buildGSoln function
    lsq = l*l
    etasq = grid.rsq/lsq
    psi_r = (-1/(2*np.pi*grid.rsq) *
             ((1 - np.exp(-etasq*etasq)) +
              2*etasq*etasq*np.log(grid.rsq)*np.exp(-etasq*etasq)))
    return psi_r

def solvePoisson(f, grid, grow=True):
    # Function for solving (- Delta) psi = f by spectral Hermite collocation
    # f can either be a callable function or a numpy array with values
    #     f[i,j] = f(x[i], y[j])
    # WARNING: numpy.meshgrid uses the convention f[i,j] = f(x[j], y[i]) which
    # is the TRANSPOSE of the setup here
    M = grid.xvect().size
    N = grid.yvect().size
    if hasattr(f, '__call__'):
        rhs = f(grid.xmat,grid.ymat)
    else:
        rhs = f

    L = -1.0 * build2dPoissonMat(grid.xvect(), grid.yvect())

    if grow == True:
        psiG, zetaG = buildGSoln(grid)
        rhs = rhs - zetaG
        psiD = np.linalg.solve(L, rhs.reshape(M*N,)).reshape(M,N)
        psi = psiD + psiG
    else:
        psi = np.linalg.solve(L, rhs.reshape((M*N,))).reshape((M,N))
    return psi

def solveScreenedPoisson(f, grid, defrad):
    # Function for solving
    #    (-Delta) psi + 1/defrad**2 psi = f
    # by spectral Hermite collocation.  f can either be a callable
    # function or a numpy array with values
    #     f[i,j] = f(x[i], y[j])
    # WARNING: numpy.meshgrid uses the convention f[i,j] = f(x[j], y[i]) which
    # is the TRANSPOSE of the setup here
    M = grid.xvect().size
    N = grid.yvect().size

    if hasattr(f, '__call__'):
        rhs = f(grid.xmat, grid.ymat)
    else:
        rhs = f

    L = build2dScreenedPoissonMat(grid.xvect(), grid.yvect(), defrad)
    psi = np.linalg.solve(L, rhs.reshape((M*N,))).reshape((M,N))
    return psi

def solveTwoLayerPoissonDirect(q1, q2, grid, defrad, grow=True):
    # Function for solving the system
    #    (-Delta) psi1 + 1/defrad**2(psi1 - psi2) = f1
    #    (-Delta) psi2 - 1/defrad**2(psi1 - psi2) = f1
    # by spectral Hermite collocation.  f can either be a callable
    # function or a numpy array with values
    #     f[i,j] = f(x[i], y[j])
    # WARNING: numpy.meshgrid uses the convention f[i,j] = f(x[j], y[i]) which
    # is the TRANSPOSE of the setup here
    M = grid.xvect().size
    N = grid.yvect().size

    if hasattr(q1, '__call__'):
        r1 = q1(grid.xmat, grid.ymat)
    else:
        r1 = q1

    if hasattr(q2, '__call__'):
        r2 = q2(grid.xmat, grid.ymat)
    else:
        r2 = q2
    rhs = np.zeros(2*M*N)
    rhs[:M*N] = r1.reshape(r1.size)
    rhs[M*N:] = r2.reshape(r2.size)

    L = build2dTwoLayerPoissonMat(grid.xvect(), grid.yvect(), defrad)
    if grow == True:
        psiG, zetaG = buildGSoln(grid)
        if not isinstance(grid, TwoLayerGrid):
            zetaG = np.tile(zetaG, (2,1,1))
            psiG = np.tile(psiG, (2,1,1))
        rhs = rhs - zetaG.reshape(zetaG.size)
        psiD = np.linalg.solve(L, rhs).reshape((2,M,N))
        psi = psiD + psiG
    else:
        psi = np.linalg.solve(L, rhs).reshape((2,M,N))
    return psi[0], psi[1]

def solveTwoLayerPoisson(q1, q2, grid, defrad, grow=True):
    # Function for solving the system
    #    (-Delta) psi1 + 1/defrad**2(psi1 - psi2) = f1
    #    (-Delta) psi2 - 1/defrad**2(psi1 - psi2) = f1
    # by spectral Hermite collocation.  f can either be a callable
    # function or a numpy array with values
    #     f[i,j] = f(x[i], y[j])
    # WARNING: numpy.meshgrid uses the convention f[i,j] = f(x[j], y[i]) which
    # is the TRANSPOSE of the setup here
    q_bt = 0.5*(q1 + q2)
    q_bc = 0.5*(q1 - q2)
    psi_bt = solvePoisson(q_bt, grid, grow=grow)
    psi_bc = solveScreenedPoisson(q_bc, grid, np.sqrt(0.5) * defrad)
    psi1 = psi_bt + psi_bc
    psi2 = psi_bt - psi_bc
    return psi1, psi2

def solveTwoLayerMatrices(q1, q2, grid, Dbt, Dbc, grow=True):
    # Function for solving the system
    #    (-Delta) psi1 + 1/defrad**2(psi1 - psi2) = f1
    #    (-Delta) psi2 - 1/defrad**2(psi1 - psi2) = f1
    # by spectral Hermite collocation.  f can either be a callable
    # function or a numpy array with values
    #     f[i,j] = f(x[i], y[j])
    # WARNING: numpy.meshgrid uses the convention f[i,j] = f(x[j], y[i]) which
    # is the TRANSPOSE of the setup here
    M,N = grid.xmat.shape

    q_bt = 0.5*(q1 + q2)
    q_bc = 0.5*(q1 - q2)

    if grow == True:
        psiG, zetaG = buildGSoln(grid)
        rhs = q_bt - zetaG
        psiD = np.linalg.solve(Dbt, rhs.reshape(M*N,)).reshape(M,N)
        psi_bt = psiD + psiG
    else:
        psi_bt = np.linalg.solve(Dbt, q_bt.reshape((M*N,))).reshape((M,N))

    psi_bc = np.linalg.solve(Dbc, q_bc.reshape((M*N,))).reshape((M,N))

    psi1 = psi_bt + psi_bc
    psi2 = psi_bt - psi_bc
    return psi1, psi2
