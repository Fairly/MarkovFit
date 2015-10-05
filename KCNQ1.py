import math

__author__ = 'fairly'

import numpy as np
import scipy.integrate as si
import matplotlib.pyplot as plt


def get_v(t):
    if type(t) == np.ndarray:
        v = [get_v(it) for it in t]
        v = np.array(v)
    else:
        if t < 3:
            v = -80
        elif 3 <= t <= 7:
            v = 20
        else:
            v = -40
    return v


F = 96485.3365
R = 8314
T = 310
RTF = R*T/F

alpha = 9.57E-4 * np.exp(1.98E-1 - Vm * RTF)
beta = 5.00E-5 * np.exp(-3.33E-2 * Vm * RTF)
gama = 3.77E-2 * np.exp(3.33E-2 * Vm * RTF)
delta = 4.77E-2 * np.exp(-4.77E-1 * Vm * RTF)
eta = 2.06E-2 * np.exp(-1.04 * Vm * RTF)
theta = 7.98E-2
psi = 2.54E-2 * np.exp(6.46E-2 * Vm * RTF)
omega = 1.78E-2 * np.exp(-5.32E-1 * Vm * RTF)
lamda = 2.32E-2 * np.exp(1.21E-1 * Vm * RTF)
mu = 6.19E-2 * np.exp(-9.46E-2 * Vm * RTF)
v = -9.00E1
x_psi = np.exp(z_psi * (z_psi+z_omega) * v*RTF)
x_omega = np.exp(z_omega * (z_psi+z_omega) * v*RTF)
psi_n = psi * x_psi**(n-1)
omega_n = omega * x_omega**(n-1)

A = np.array([
    [0, 3.8, 0],
    [-4.9, 0, 1.0],
    [0, -0.67, 0]
])
B = np.array([
    [0, 0.028, 0],
    [0.008, 0, 0.041],
    [0, -0.06, 0]
])


def f(y, t0, A, B):
    c0, c1, o2 = y

    V = get_v(t0)
    k01 = np.exp(A[0][1] + B[0][1] * V)
    k10 = np.exp(A[1][0] + B[1][0] * V)
    k12 = np.exp(A[1][2] + B[1][2] * V)
    k21 = np.exp(A[2][1] + B[2][1] * V)

    dc0 = k10 * c1 - k01 * c0
    dc1 = k01 * c0 + k21 * o2 - (k10 + k12) * c1
    do2 = k12 * c1 - k21 * o2

    return [dc0, dc1, do2]


def test_KvLQT1():
    t = np.linspace(0, 10.0, num=10001)
    y0 = [1, 0, 0]
    y = si.odeint(f, y0, t, args=(A, B))

    G_ = 12
    R = 8314.472
    T = 295.15
    z = 1
    F = 96485.3365
    Ko = 4.0
    Ki = 140
    RTzF = (R * T)/(z * F)
    Ek = RTzF * math.log(Ko/Ki)
    I = G_ * y[:, 2] * (get_v(t) - Ek)
    print I

    plt.figure(1)
    plt.plot(t, I, '--', label='Ik')
    plt.legend()

    plt.figure(2)
    plt.plot(t, y[:, 0], 'r', label='C0')
    plt.plot(t, y[:, 1], 'g', label='C1')
    plt.plot(t, y[:, 2], 'b', label='O2')
    plt.axis([1.0, 8, -0.1, 1.0])
    plt.legend()
    plt.show()


def print_k():
    for V in range(-80, 60, 10):
        k01 = np.exp(A[0][1] + B[0][1] * V)
        k10 = np.exp(A[1][0] + B[1][0] * V)
        k12 = np.exp(A[1][2] + B[1][2] * V)
        k21 = np.exp(A[2][1] + B[2][1] * V)
        print(V, [k01, k10, k12, k21])


if __name__ == '__main__':
    # print_k()
    test_KvLQT1()
    pass