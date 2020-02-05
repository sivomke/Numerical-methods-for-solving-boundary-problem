# implementing Bubnov-Galerkin method for solving boundary value problem with Robin conditions
import diff_system as sys  # system definition
import coord_system as coord  # definition of basis functions
import adaptive as ad
import numpy as np
import scipy
import plotting as plt

n = 12  # number of basis functions
tol = 1e-5


def approx_decorator(coef: np.ndarray):
    return lambda x: sum(coef[i - 1]*coord.phi(i, x) for i in range(1, n+1)) + sys.A*x + sys.B


def bubnov_galerkin():
    G = np.array([[ad.adaptive(lambda x: coord.A_phi(i, x)*coord.phi(j, x), sys.a, sys.b, tol) for i in range(1, n + 1)]
                  for j in range(1, n + 1)])
    y = np.array([ad.adaptive(lambda x: sys.f1(x) * coord.phi(i, x), sys.a, sys.b, tol) for i in range(1, n + 1)])
    # G = np.array([[scipy.integrate.quad(lambda x: A_phi(i, x)*phi(j, x), a, b)[0] for i in range(1, n + 1)] for j in
    #              range(1, n + 1)])
    # y = np.array([scipy.integrate.quad(lambda x: f1(x) * phi(i, x), a, b)[0] for i in range(1, n + 1)])
    coef = np.array(np.linalg.solve(G, y))
    print("Residual coefficients: ", np.dot(G, coef) - y)
    plt.comparison(approx_decorator(coef), sys.u, sys.a, sys.b, -10, 10)


# bubnov_galerkin()