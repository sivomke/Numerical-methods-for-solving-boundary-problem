import numpy as np
import math


# system parameters
# -(k(x)*u'(x))' + p(x)*u'(x) + q(x)*u(x) = f(x), a<x<b
# -k(x)*u'(x) + alpha_1*u(x) = mu_1
# k(x)*u'(x) + alpha_2*u(x) = mu_2

a = 0
b = 3
m1 = 3
m2 = 5
m3 = 2
m4 = 2
m5 = 1
p1 = 2
p2 = 2
p3 = 1
q1 = 0
q2 = 2
q3 = 1
k1 = 2
k2 = 2
k3 = 1
alpha_1 = 1
alpha_2 = 3


# exact system solution
def u(x):
    return m1*math.sin(m2*x) + m3


def k(x):
    return k1*math.sin(k2*x) + k3


def p(x):
    return p1*math.sin(p2*x) + p3


def q(x):
    return q1*math.sin(q2*x) + q3


def dk_dx(x):
    return k1*k2*math.cos(k2*x)


def du_dx(x):
    return m1*m2*math.cos(m2*x)


def d2u_dx2(x):
    return -m1*m2**2*math.sin(m2*x)


def f(x):
    return -(dk_dx(x)*du_dx(x)+k(x)*d2u_dx2(x)) + p(x)*du_dx(x) + q(x)*u(x)


# boundary conditions
mu_1 = -k(a)*du_dx(a) + alpha_1*u(a)  # (1)
mu_2 = k(b)*du_dx(b) + alpha_2*u(b)  # (2)


# we want to get system with homogeneous boundary conditions, so we present u(x) = v(x) + psi(x), where
# psi(x) satisfies inhomogeneous boundary conditions
# psi(x) = A*x + B, where const A,B we find solving (1),(2) with respect to A and B

M = np.array([[alpha_1*a - k(a), alpha_1], [alpha_2*b + k(b), alpha_2]])
mu_vec = np.array([mu_1, mu_2])
A, B = np.linalg.solve(M, mu_vec)
# print("A: ", A, "B: ", B)
# print("Residual for A and B: ", np.dot(M, np.array([A, B])) - mu_vec)


# we obtain new problem with respect to v(x) with new right hand-side f1(x) and homogeneous boundary conditions
def f1(x):
    return f(x) + dk_dx(x)*A - p(x)*A - q(x)*(A*x + B)