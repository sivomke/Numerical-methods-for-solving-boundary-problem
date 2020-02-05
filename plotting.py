import numpy as np
import matplotlib.pyplot as plt


# plotting target func and its numeric approx
def comparison(approx, target, left, right, bottom, top):
    t = np.linspace(left, right, 1000)
    approx_val = [approx(time) for time in t]
    target_val = [target(time) for time in t]
    plt.plot(t, approx_val, label='approximation function')
    plt.plot(t, target_val, label='target function')
    plt.axis([left, right, bottom, top])
    plt.legend()
    plt.show()


def plot(func, left, right, bottom, top):
    t = np.linspace(left, right, 1000)
    func_val = [func(time) for time in t]
    plt.plot(t, func_val)
    plt.axis([left, right, bottom, top])
    plt.xlabel("frequency, Hz")
    plt.show()