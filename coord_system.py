import diff_system as sys  # definition of diff equations system


# system of functions phi(1,x) = (x-sys.a)**2*(x-C), phi(2,x) = (b-x)**2*(x-D), phi(i,x) = (x-a)**(i-1)*(b-x)**2, i=3..n
# in order to this system of functions to be complete we choose C such that phi(1,x) satisfies (2)
# and we choose D such that phi(2,x) satisfies (2)
C = sys.b+sys.k(sys.b)*(sys.b-sys.a)/(2*sys.k(sys.b)+sys.alpha_2*(sys.b-sys.a))
D = sys.a-sys.k(sys.a)*(sys.b-sys.a)/(2*sys.k(sys.a)+sys.alpha_1*(sys.b-sys.a))
print("C: ", C)
print("D: ", D)


# system of functions
def phi(i, x):
    if i == 1:
        return (x-sys.a)**2*(x-C)
    elif i == 2:
        return (sys.b-x)**2*(x-D)
    else:
        return (x-sys.a)**(i-1)*(sys.b-x)**2


def dphi_dx(i, x):
    if i == 1:
        return 2*(x-sys.a)*(x-C) + (x-sys.a)**2
    elif i == 2:
        return -2*(sys.b-x)*(x-D) + (sys.b-x)**2
    else:
        return (i-1)*(x-sys.a)**(i-2)*(sys.b-x)**2-2*(sys.b-x)*(x-sys.a)**(i-1)


def d2phi_dx2(i, x):
    if i == 1:
        return 2*(x-C) + 4*(x-sys.a)
    elif i == 2:
        return 2*(x-D) - 4*(sys.b-x)
    else:
        return (i-1)*(i-2)*(x-sys.a)**(i-3)*(sys.b-x)**2 - 4*(sys.b-x)*(i-1)*(x-sys.a)**(i-2) + 2*(x-sys.a)**(i-1)


def A_phi(i, x):
    return dphi_dx(i, x)*(-sys.dk_dx(x) + sys.p(x)) - sys.k(x)*d2phi_dx2(i, x) + sys.q(x)*phi(i, x)