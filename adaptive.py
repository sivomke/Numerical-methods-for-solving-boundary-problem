def simpson(func, left, right):
    return (right - left)*(func(left) + 4*func(0.5*(left + right)) + func(right))/6


# here we check the tolerance by Runge rule
# the result is enhanced by Richardson extrapolation formula
# p = 4 (accuracy by step)
def adaptive_simpson(func, left, right, tol):
    mid = 0.5*(left + right)
    s1 = simpson(func, left, right)
    s2 = simpson(func, left, mid) + simpson(func, mid, right)
    if abs(s2 - s1) < 15*tol:
        return s2 + (s2 - s1)/15
    else:
        return adaptive_simpson(func, left, mid, 0.5*tol) + adaptive_simpson(func, mid, right, 0.5*tol)


# here we divide our main interval into N sub-intervals and use adaptive Simpson method on each of those sub-intervals
def adaptive(func, left, right, tol):
    num = 200
    step = (right - left)/num
    return sum(adaptive_simpson(func, left + (i-1)*step, left + i*step, tol/num) for i in range(1, num + 1))
