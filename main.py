import math
import numpy as np
import xlrd
import scipy.integrate as inte


def xls():
    data = xlrd.open_workbook('steady-state-I-V.xlsx')
    table = data.sheet_by_index(1)
    c = table.ncols
    r = table.nrows
    data_points = 1
    print(np.array(table.col_values(0), table.col_values(1)))


def get_v(t):
    if t < 0:
        v = -80
    elif 0 <= t <= 2000:
        v = 20
    else:
        v = -40
    return v


def model(x, k):
    # physical constants
    R = 8314.0
    T = 310.0
    F = 96485.0
    # inte.odeint(func, y0=, t)



if __name__ == "__main__":
    # initialize
    k = np.zeros((3, 3))
    a = np.zeros((3, 3))
    b = np.zeros((3, 3))
    v = 0
    t = 0

    a[0, 1] = 3.8
    a[1, 0] = -4.9
    a[1, 2] = 1.0
    a[2, 1] = -0.67

    b[0, 1] = 0.028
    b[1, 0] = 0.008
    b[1, 2] = 0.041
    b[2, 1] = -0.06

    for i in xrange(3):
        for j in xrange(3):
            k[i, j] = math.exp(a[i, j] + b[i, j] * v)

