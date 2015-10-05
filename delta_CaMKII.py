__author__ = 'fairly'

import pprint

import xlrd

import numpy as np
from scipy.integrate import odeint
from scipy.optimize import leastsq
from scipy.optimize import minimize
import matplotlib.pyplot as plt


def get_Cai(t0):
    t = t0 % (1000 / 3)
    if 30 <= t < 90:
        return 2.5E-5 * (t - 30)
    elif 90 <= t < 200:
        return 1.5E-3 * (200 - t) / 110
    else:
        return 0


def CaM_transition(y, t0, CaM_total):
    Cai = get_Cai(t0)

    # CaM_total = 0.006
    CaMCa, CaMCa2, CaMCa3, CaMCa4 = y
    CaM = CaM_total - CaMCa - CaMCa2 - CaMCa3 - CaMCa4

    k_CaCaM_1 = 2.5
    k_CaCaM_1_minus = 0.05
    k_CaCaM_2 = 88.25
    k_CaCaM_2_minus = 0.05
    k_CaCaM_3 = 12.5
    k_CaCaM_3_minus = 1.25
    k_CaCaM_4 = 250
    k_CaCaM_4_minus = 1.25

    CA1 = -k_CaCaM_1 * Cai * CaM + k_CaCaM_1_minus * CaMCa
    CA2 = -k_CaCaM_2 * Cai * CaMCa + k_CaCaM_2_minus * CaMCa2
    CA3 = -k_CaCaM_3 * Cai * CaMCa2 + k_CaCaM_3_minus * CaMCa3
    CA4 = -k_CaCaM_4 * Cai * CaMCa3 + k_CaCaM_4_minus * CaMCa4

    d1 = -CA1 + CA2
    d2 = -CA2 + CA3
    d3 = -CA3 + CA4
    d4 = -CA4

    return [d1, d2, d3, d4]


def CaMKII_transition(y, t0, CaM_total, CaMKII_total, Cai, ATP, temperature,
                      param, frequency=False):
    """

    :param y:
    :param t0:
    :param CaM_total:
    :param CaMKII_total:
    :param Cai:
    :param ATP:
    :param temperature:
    :param param:
    :param frequency: This parameter is designed for reproducing the Figure4, whose protocol is different with others.
    :return:
    """
    CaMKII, CaMKII_CaMCa4, CaMKIIP_CaMCa4, CaMCa, CaMCa2, CaMCa3, CaMCa4 = y
    CaMKIIP = CaMKII_total - CaMKII - CaMKII_CaMCa4 - CaMKIIP_CaMCa4
    # TODO original equation modified below
    CaM = CaM_total - CaMCa - CaMCa2 - CaMCa3 - CaMCa4 - CaMKII_CaMCa4 - CaMKIIP_CaMCa4
    # CaM = CaM_total - CaMCa - CaMCa2 - CaMCa3 - CaMCa4

    if not frequency:
        pass
    else:
        # exposed in mixture fluid for a 200ms duration in each cycle
        if t0 % np.floor(1000 / frequency) <= 200:
            pass
        else:
            Cai = 0

    k_CaCaM_1 = 2.5
    k_CaCaM_1_minus = 0.05
    k_CaCaM_2 = 88.25
    k_CaCaM_2_minus = 0.05
    k_CaCaM_3 = 12.5
    k_CaCaM_3_minus = 1.25
    k_CaCaM_4 = 250
    k_CaCaM_4_minus = 1.25

    CA1 = -k_CaCaM_1 * Cai * CaM + k_CaCaM_1_minus * CaMCa
    CA2 = -k_CaCaM_2 * Cai * CaMCa + k_CaCaM_2_minus * CaMCa2
    CA3 = -k_CaCaM_3 * Cai * CaMCa2 + k_CaCaM_3_minus * CaMCa3
    CA4 = -k_CaCaM_4 * Cai * CaMCa3 + k_CaCaM_4_minus * CaMCa4

    kcat, k_asso, k_disso, k_dissoCa = param
    # kcat = args
    # k_asso = 2.1
    # k_disso = 0.7E-4
    # k_dissoCa = 0.95E-3
    k_disso2 = k_disso / 1000
    k_dissoCa2 = k_dissoCa / 1000
    KmCaM = 3.0E-5

    if temperature == 30:
        kcat *= 30
    if temperature == 37:
        kcat *= 90

    KmATP = 19.1E-3
    kcat_PP1 = 1.72E-3
    Km_PP1 = 11.0E-3
    PP1 = 0  # the concentration of PP1 in figure 3 is 0

    A1 = k_asso * CaMCa4 * CaMKII
    A2 = (k_disso * (1 - KmCaM ** 3 / (Cai ** 3 + KmCaM ** 3)) + k_dissoCa * (
        KmCaM ** 3 / (Cai ** 3 + KmCaM ** 3))) * CaMKII_CaMCa4
    P = 1 - (CaMKII / CaMKII_total) ** 2
    B1 = kcat * P * (ATP / (ATP + KmATP)) * CaMKII_CaMCa4
    B2 = kcat_PP1 * CaMKIIP_CaMCa4 / (CaMKIIP_CaMCa4 + Km_PP1) * PP1
    D1 = kcat_PP1 * CaMKIIP / (CaMKIIP_CaMCa4 + Km_PP1) * PP1
    C1 = (k_disso2 * (1 - KmCaM ** 3 / (Cai ** 3 + KmCaM ** 3)) + k_dissoCa2 * (
        KmCaM ** 3 / (Cai ** 3 + KmCaM ** 3))) * CaMKIIP_CaMCa4
    C2 = k_asso * CaMKIIP * CaMCa4

    dCaMKII = -A1 + A2 + D1
    dCaMKII_CaMCa4 = A1 - A2 - B1 + B2
    dCaMKIIP_CaMCa4 = B1 - B2 - C1 + C2
    d1 = -CA1 + CA2
    d2 = -CA2 + CA3
    d3 = -CA3 + CA4
    d4 = -CA4 - A1 + k_disso * (1 - KmCaM ** 3 / (Cai ** 3 + KmCaM ** 3)) * CaMKII_CaMCa4 \
         - C2 + k_disso2 * (1 - KmCaM ** 3 / (Cai ** 3 + KmCaM ** 3)) * CaMKIIP_CaMCa4  # TODO
    # d4 = -CA4

    return [dCaMKII, dCaMKII_CaMCa4, dCaMKIIP_CaMCa4, d1, d2, d3, d4]


def CaMKIIROS_transition(y, t0, CaM_total, CaMKII_total, Cai, ATP, temperature, ROS, MsrA,
                         param, frequency=False):
    """

    :param y:
    :param t0:
    :param CaM_total:
    :param CaMKII_total:
    :param Cai:
    :param ATP:
    :param temperature:
    :param param:
    :param frequency: This parameter is designed for reproducing the Figure4, whose protocol is different with others.
    :return:
    """
    CaMKII, CaMKII_CaMCa4, CaMKIIP_CaMCa4, OX, OXA, CaMCa, CaMCa2, CaMCa3, CaMCa4 = y
    CaMKIIP = CaMKII_total - CaMKII - CaMKII_CaMCa4 - CaMKIIP_CaMCa4 - OX - OXA
    CaM = CaM_total - CaMCa - CaMCa2 - CaMCa3 - CaMCa4 - CaMKII_CaMCa4 - CaMKIIP_CaMCa4 - OX

    if not frequency:
        pass
    else:
        # exposed in mixture fluid for a 200ms duration in each cycle
        if t0 % np.floor(1000 / frequency) <= 200:
            pass
        else:
            Cai = 0

    k_CaCaM_1 = 2.5
    k_CaCaM_1_minus = 0.05
    k_CaCaM_2 = 88.25
    k_CaCaM_2_minus = 0.05
    k_CaCaM_3 = 12.5
    k_CaCaM_3_minus = 1.25
    k_CaCaM_4 = 250
    k_CaCaM_4_minus = 1.25

    CA1 = -k_CaCaM_1 * Cai * CaM + k_CaCaM_1_minus * CaMCa
    CA2 = -k_CaCaM_2 * Cai * CaMCa + k_CaCaM_2_minus * CaMCa2
    CA3 = -k_CaCaM_3 * Cai * CaMCa2 + k_CaCaM_3_minus * CaMCa3
    CA4 = -k_CaCaM_4 * Cai * CaMCa3 + k_CaCaM_4_minus * CaMCa4

    kcat, k_asso, k_disso, k_dissoCa = [1.99161904e-04, 1.00652800e+00, 2.46811261e-06,
                                        1.27324832e-03]
    k_ox, = param
    # kcat = args
    # k_asso = 2.1
    # k_disso = 0.7E-4
    # k_dissoCa = 0.95E-3
    k_disso2 = k_disso / 1000
    k_dissoCa2 = k_dissoCa / 1000
    # k_ox = 1.952869E-6
    k_red = 0.28   # this one and the Km_MsrA is from literature

    if temperature == 30:
        kcat *= 30
    if temperature == 37:
        kcat *= 90

    KmCaM = 3.0E-5
    KmATP = 19.1E-3
    kcat_PP1 = 1.72E-3
    Km_PP1 = 11.0E-3
    PP1 = 0  # the concentration of PP1 in figure 3 is 0
    Km_ROS = 0.06E-3
    Km_MsrA = 0.34

    cons_cam = (1 - KmCaM ** 3 / (Cai ** 3 + KmCaM ** 3))

    A1 = k_asso * CaMCa4 * CaMKII
    A2 = (k_disso * cons_cam + k_dissoCa * (KmCaM ** 3 / (Cai ** 3 + KmCaM ** 3))) * CaMKII_CaMCa4
    P = 1 - (CaMKII / CaMKII_total) ** 2
    B1 = kcat * P * (ATP / (ATP + KmATP)) * CaMKII_CaMCa4
    B2 = kcat_PP1 * CaMKIIP_CaMCa4 / (CaMKIIP_CaMCa4 + Km_PP1) * PP1
    D1 = kcat_PP1 * CaMKIIP / (CaMKIIP_CaMCa4 + Km_PP1) * PP1
    C1 = (k_disso2 * cons_cam + k_dissoCa2 * (KmCaM ** 3 / (Cai ** 3 + KmCaM ** 3))) * CaMKIIP_CaMCa4
    C2 = k_asso * CaMKIIP * CaMCa4
    E1 = k_ox * ROS / (ROS + Km_ROS) * CaMKII_CaMCa4
    E2 = k_red * MsrA / (MsrA + Km_MsrA) * OX
    G2 = k_red * MsrA / (MsrA + Km_MsrA) * OXA
    F1 = (k_disso * cons_cam + k_dissoCa * (KmCaM ** 3 / (Cai ** 3 + KmCaM ** 3))) * OX
    F2 = k_asso * OXA * CaMCa4

    dCaMKII = -A1 + A2 + D1 + G2
    dCaMKII_CaMCa4 = A1 - A2 - B1 + B2 + E2 - E1
    dCaMKIIP_CaMCa4 = B1 - B2 - C1 + C2
    dOX = E1 + F2 - E2 - F1
    dOXA = -G2 - F2 + F1
    d1 = -CA1 + CA2
    d2 = -CA2 + CA3
    d3 = -CA3 + CA4
    d4 = -CA4 - A1 + k_disso * cons_cam * CaMKII_CaMCa4 - C2 + k_disso2 \
                                                               * cons_cam * CaMKIIP_CaMCa4 + k_disso * cons_cam * OX \
         - F2
    # d4 = -CA4

    return [dCaMKII, dCaMKII_CaMCa4, dCaMKIIP_CaMCa4, dOX, dOXA, d1, d2, d3, d4]


def curvfit(data, nomalized=True):
    a = 1
    b = 1
    h = 1
    param = [a, b, h]

    fit_param, flag = leastsq(err_equation, param, args=(data,))

    if nomalized:
        # normalized to the maximum of the curve, just for the figures having a percentage y axis
        maximum = fit_hillequation(data[-1][0], *fit_param)
        for i, item in enumerate(data):
            data[i][1] /= maximum

    fit_param, flag = leastsq(err_equation, param, args=(data,))

    return fit_param


def err_equation(param, data):
    err = data[:, 1] - fit_hillequation(data[:, 0], *param)
    # print sum(err**2)
    return err


def err_figure2(param, data_2a, data_2b, data_2c):
    t0 = np.linspace(0, 60000, num=60001)  # all 3 simulations last 1 minute
    ATP = 0

    CaMKII_total = 5E-6
    CaMKII_initial = CaMKII_total
    Cai = 0.5
    y0 = [CaMKII_initial, 0, 0, 0, 0]
    # error_2a = np.zeros(len(data_2a))
    #
    # for i, (CaM_total, fraction) in enumerate(data_2a):  # calculate the expected data points
    # y = odeint(CaMKII_transition_noatp, y0, t0, args=(CaM_total, CaMKII_total, Cai, param))
    # error_2a[i] = fraction-(CaMKII_total-y[-1][0])/CaMKII_total

    error_2b = np.zeros(len(data_2b))
    CaM_total = 1E-4  # CaM_total is set to 1E-4mM
    Cai = 0.5
    for i, (CaMKII_total, fraction) in enumerate(data_2b):
        y0 = [CaMKII_total, 0, 0, 0, 0, 0, 0]
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, param))
        error_2b[i] = fraction - (CaMKII_total - y[-1][0]) / CaM_total

    # error_2c = np.zeros(len(data_2c))
    # CaMKII_total = 2E-4
    # CaM_total = 5E-3
    # y0 = [CaMKII_total, 0, 0, 0, 0]
    # for i, (Cai, fraction) in enumerate(data_2c):
    # y = odeint(CaMKII_transition_noatp, y0, t0, args=(CaM_total, CaMKII_total, Cai, param))
    # error_2c[i] = fraction - (CaMKII_total-y[-1][0])/CaMKII_total

    # error = np.hstack((error_2a, error_2b, error_2c))
    error = error_2b
    print sum(error ** 2)
    return error


def err_figure3(param, data_3a, data_3b, data_3c):
    # t0 = np.linspace(0, 15000, num=15001)
    # CaMKII_total = 62E-6
    # CaMKII_initial = CaMKII_total
    # Cai = 0.5
    # ATP = 0.1
    # y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0]
    # error_3a = np.zeros(len(data_3a))
    #
    # for i, (CaM_total, percentage) in enumerate(data_3a):  # calculate the expected data points
    # y = odeint(CaMKII_transition_atp, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, param))
    # error_3a[i] = percentage-100*(CaMKII_total-y[-1][0]-y[-1][1])/CaMKII_total

    t0 = np.linspace(0, 6000, num=6001)
    CaMKII_total = 5E-6
    CaMKII_initial = CaMKII_total
    Cai = 0.5
    ATP = 0.25
    y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0]
    error_3b = np.zeros(len(data_3b))

    for i, (CaM_total, percentage) in enumerate(data_3b):
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, param))
        error_3b[i] = percentage - 100 * (CaMKII_total - y[-1][0] - y[-1][1]) / CaMKII_total

    t0 = np.linspace(0, 300000, num=300001)
    CaMKII_total = 0.2E-3
    CaMKII_initial = CaMKII_total
    CaM_total = 50E-3
    ATP = 2
    y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0]
    error_3c = np.zeros(len(data_3c))

    for i, (Cai, fraction) in enumerate(data_3c):
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, param))
        error_3c[i] = fraction - 100 * (CaMKII_total - y[-1][0] - y[-1][1]) / CaMKII_total

    error = np.hstack((error_3b, error_3c))  # TODO
    # error = error_3a
    print sum(error ** 2)
    return error


def err_figure4(param, curve_41, curve_42, curve_44):
    # 1Hz error
    x = np.linspace(0, 100000, num=100001)  # simulate 100s for the 1Hz case
    CaMKII_total = 5E-6  # because the CaMKII is fixed in a tube where fluid is flowing,
    # so the concentration is estimated
    CaM_total = 0.1
    Cai = 0.5
    ATP = 0.25
    y0 = [CaMKII_total, 0, 0, 0, 0, 0, 0]
    error_41 = np.zeros(len(x))

    # for i, (CaM_total, percentage) in enumerate(x):
    # y = odeint(CaMKII_transition_atp, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, param))
    # error_3b[i] = percentage - 100 * (CaMKII_total - y[-1][0] - y[-1][1]) / CaMKII_total
    # pass


def err_figure5(param, curve_5a, curve_5b, curve_5c):
    t0 = np.linspace(0, 60000, num=60001)
    ATP = 0
    Cai = 0.5
    temperature = 30
    CaM_total = 1E-4  # CaM_total is set to 1E-4mM
    x = np.logspace(-6, -2, num=30)
    error_5a = np.zeros(len(x))

    for i, CaMKII_total in enumerate(x):
        y0 = [CaMKII_total, 0, 0, 0, 0, 0, 0]
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param))
        error_5a[i] = fit_hillequation(CaMKII_total, *curve_5a) - (CaMKII_total - y[-1][0]) / CaM_total

    t0 = np.linspace(0, 15000, num=15001)
    CaMKII_total = 62E-6
    CaMKII_initial = CaMKII_total
    Cai = 0.5
    ATP = 0.1
    temperature = 0
    y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0]
    x = np.logspace(-6, -2, num=30)
    error_5b = np.zeros(len(x))

    for i, CaM_total in enumerate(x):
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param))
        error_5b[i] = (CaMKII_total - y[-1][0] - y[-1][1]) / CaMKII_total

    for i, CaM_total in enumerate(x):
        error_5b[i] = fit_hillequation(CaM_total, *curve_5b) - error_5b[i] / error_5b[-1]

    t0 = np.linspace(0, 60000, num=60001)
    CaMKII_total = 6.2E-6
    CaMKII_initial = CaMKII_total
    Cai = 0.5
    ATP = 0.1
    temperature = 30
    y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0]
    x = np.logspace(-7, -2, num=30)
    error_5c = np.zeros(len(x))

    for i, CaM_total in enumerate(x):
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param))
        error_5c[i] = fit_hillequation(CaM_total, *curve_5c) - (CaMKII_total - y[-1][0]) / CaMKII_total

    error = np.hstack((error_5a, error_5b, error_5c))
    # error = error_5a
    print sum(error ** 2)
    return sum(error ** 2)
    # return error


def err_ROS(param, data_ROS):
    t0 = np.linspace(0, 60000, num=60001)
    CaMKII_total = 6.2E-4
    CaMKII_initial = CaMKII_total
    CaM_total = 1E-3
    Cai = 200E-3
    ATP = 0
    temperature = 0
    y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0, 0, 0]
    x = [0, 0.01E-3, 0.1E-3, 1E-3]
    MsrA = 0
    result = np.zeros((len(x), 2))

    tmp_max = max([k_a for k_a in data_ROS[:, 1]])
    data_ROS = np.array([[ROS, k_a / tmp_max * 100] for ROS, k_a in data_ROS])

    for i, ROS in enumerate(x):
        y = odeint(CaMKIIROS_transition, y0, t0, args=(
            CaM_total, CaMKII_total, Cai, ATP, temperature, ROS, MsrA, param))
        # result[i] = [ROS, (CaMKII_total - y[-1][0]) / CaMKII_total]
        # tmp = [y[-1][0], y[-1][1], y[-1][2],
        #        CaMKII_total - y[-1][0] - y[-1][1] - y[-1][2] - y[-1][3] - y[-1][4], y[-1][3], y[-1][4]]
        # TODO
        # notice this simplification, due to the experimental procedures,
        # all the kinase will turn to be OXA at the end, so I just count the OX and OXA is OK
        tmp = y[-1][3:5]
        result[i] = [ROS, sum(tmp)]

    result_m = max([i for a, i in result])
    # result = np.array([[ROS, k_a / result[-1][1]] for ROS, k_a in result])
    result = np.array([[ROS, k_a / result_m * 100] for ROS, k_a in result])

    err = data_ROS[:, 1] - result[:, 1]
    err = sum(err ** 2)
    print err
    return err


def fit_figure(*figures):
    data = [[name, np.array(xls(name))] for name in figures]

    param_alpha = np.array([2.1, 1.4E-4, 1.9E-3, 1E-5])  # todo 1000 times
    param_delta = np.array([2.1, 0.7E-4, 0.95E-3, 6E-5])

    for fit in data:
        if fit[0].startswith('5'):
            fit.append(param_delta)
        else:
            fit.append(param_alpha)

    data = set_environment(data)
    print data


def fit_figure2():
    data_2a = np.array(xls('2a'))
    data_2a[:, 0] *= 1E-6  # unify the dimension to mM-1ms-1
    data_2b = np.array(xls('2b'))
    data_2b[:, 0] *= 1E-6
    data_2c = np.array(xls('2c'))
    data_2c[:, 0] *= 1E-3

    # param = np.array([2.1, 1.4E-4, 1.9E-3])
    param = np.array([4.2, 1.4E-4, 1.9E-3, 1.4E-7, 1.9E-6, 1E-5])
    # plot_figure2a(param, data_2a)
    # plot_figure2b(param, data_2b)
    # plot_figure2c(param, data_2c)

    fit_param = leastsq(err_figure2, param, args=(data_2a, data_2b, data_2c))

    print fit_param
    # plot_figure2a(fit_param[0], data_2a)
    plot_figure2b(fit_param[0], data_2b)
    # plot_figure2c(fit_param[0], data_2c)
    plt.show()


def fit_figure3():
    data_3a = np.array(xls('3a'))
    curve_3a = curvfit(data_3a)
    data_3b = np.array(xls('3b'))
    curve_3b = curvfit(data_3b)
    data_3c = np.array(xls('3c'))
    curve_3c = curvfit(data_3c)

    param = np.array([6E-5, 2.1, 1.4E-4, 1.9E-3])
    plot_figure3a(param, data_3a, curve_3a)
    # plot_figure3b(param, data_3b, curve_3b)
    # plot_figure3c(param, data_3c, curve_3c)

    # fit_param = leastsq(err_figure3, param, args=(data_3a, data_3b, data_3c))

    # print fit_param
    # plot_figure3a(fit_param[0], data_3a)
    # plot_figure3b(fit_param[0], data_3b)
    # plot_figure3c(fit_param[0], data_3c)
    plt.show()


def fit_figure4():
    data_41 = np.array(xls('1Hz'))
    curve_41 = curvfit(data_41[1:], nomalized=False)
    data_42 = np.array(xls('2.5Hz'))
    curve_42 = curvfit(data_42[1:], nomalized=False)
    data_44 = np.array(xls('4Hz'))
    curve_44 = curvfit(data_44[1:], nomalized=False)

    # param = np.array([2.73e-05, 1.15e+00, 1e-08, 2e-03])
    # param = np.array([1.49e-04, 9.8e-01, 1.89e-06, 1.09e-03])
    # param = np.array([6E-5, 2.1, 0.7E-4, 0.95E-3])
    param_a = np.array([1E-5, 2.1, 1.4E-4, 1.9E-3])
    param_d = np.array([1.49e-04, 9.8e-01, 1.89e-06, 1.09e-03])

    plot_figure4(param_a, param_d, data_41, curve_41, data_42, curve_42, data_44, curve_44)

    plt.show()


def fit_figure5():
    data_5a = np.array(xls('5a'))
    curve_5a = curvfit(data_5a)
    data_5b = np.array(xls('5b'))
    curve_5b = curvfit(data_5b)
    data_5c = np.array(xls('5c'))
    curve_5c = curvfit(data_5c)

    # param = [kcat, k_asso, k_disso, k_dissoCa]
    # alpha
    # param = np.array([1E-5, 2.1, 1.4E-4, 1.9E-3])
    # delta
    # param = np.array([6E-5, 2.1, 0.7E-4, 0.95E-3])

    # param = np.array([2.73e-05, 1.15e+00, 1e-08, 4.50e-03])
    param = np.array([1.49e-04, 9.8e-01, 1.89e-06, 1.09e-03])
    # param = np.array([2.73e-05, 1.15e+00, 1e-08, 2e-03])
    plot_figure5a(param, data_5a, curve_5a)
    plot_figure5b(param, data_5b, curve_5b)
    plot_figure5c(param, data_5c, curve_5c)

    # if need to fit, uncomment following 7 lines and comment above 3 lines.
    # fit_param = minimize(err_figure5, param, args=(curve_5a, curve_5b, curve_5c),
    #                      method='Nelder-Mead', options={'maxfev': 30})

    # print fit_param
    # plot_figure5a(fit_param['x'], data_5a, curve_5a)
    # plot_figure5b(fit_param['x'], data_5b, curve_5b)
    # plot_figure5c(fit_param['x'], data_5c, curve_5c)
    # plt.show()


def fit_ROS():
    data_ROS = xls('ROS2')

    # param = np.array([6E-5, 2.1, 0.7E-4, 0.95E-3, 1.0E-3])
    # param = np.array([2.73e-05, 1.15e+00, 1e-08, 2e-03, 1.46735612e-05])
    # param = np.array([1.952869E-6])
    param = np.array([6.48413535e-06])

    # param = np.array([2.73e-05, 1.15e+00, 1e-08, 4.50e-03])
    # param = np.array([4.24526611e-05, 5.00503650e+00, 4.83683594e-05, 9.66395965e-04])
    # param = np.array([2.73e-05, 1.15e+00, 1e-08, 2e-03])

    # fit_param = minimize(err_ROS, param, args=(data_ROS,), method='Nelder-Mead', options={'maxfev': 100})
    # param = fit_param.x

    pprint.pprint(param)
    plot_ROS(param, data_ROS)

    plt.show()


def fit_hillequation(x, a, b, h):
    return a * np.power(x, h) / (np.power(b, h) + np.power(x, h))


def get_kd(result):
    """
    linear interpolation on a logarithmic axis
    :param result:
    :return:
    """
    for i, (x, y) in enumerate(result):
        if y > 0.5:
            y1 = result[i - 1][1]
            y2 = y
            x1 = np.log10(result[i - 1][0])
            x2 = np.log10(x)
            return 10 ** ((0.5 - y1) * (x2 - x1) / (y2 - y1) + x1)


def investigation():
    data_5a = np.array(xls('5a'))
    curve_5a = curvfit(data_5a)
    data_5b = np.array(xls('5b'))
    curve_5b = curvfit(data_5b)
    data_5c = np.array(xls('5c'))
    curve_5c = curvfit(data_5c)

    x = np.linspace(0.1, 1.1, num=8)
    y = np.linspace(1E-5, 6E-3, num=8)

    z = np.zeros((len(x), len(y)))
    for i in range(0, len(x)):
        for j in range(0, len(y)):
            param = np.array([y[i], x[i], 0.7E-4, 0.95E-3])
            error = err_figure5(param, curve_5a, curve_5b, curve_5c)
            z[i][j] = sum(error ** 2)
            print i, j

    x, y = np.meshgrid(x, y)
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot_surface(x, y, z)
    plt.show()


def kinase_activity(CaMKII, CaMKII_CaMCa4=0, CaMKIIP_CaMCa4=0, CaMKIIP=0, Ox=0, OxA=0):
    args = [CaMKII, CaMKII_CaMCa4, CaMKIIP_CaMCa4, CaMKIIP, Ox, OxA]
    for i, j in enumerate(args):
        if j < sum(args) / 1000:
            args[i] = 0
    # a1 = args[0] * 5
    a2 = args[1] * 60
    a3 = (args[3] + args[2]) * 40
    a4 = (args[4] + args[5]) * 25
    # return a1 + a2 + a3 + a4
    return a2 + a3 + a4


def plot_figure2a(param, data_2a):
    t0 = np.linspace(0, 60000, num=60001)
    CaMKII_total = 5E-6
    CaMKII_initial = CaMKII_total
    Cai = 0.5
    y0 = [CaMKII_initial, 0, 0, 0, 0]
    result = np.zeros((len(data_2a), 2))

    for i, (CaM_total, fraction) in enumerate(data_2a):
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, param))
        result[i] = [CaM_total, (CaMKII_total - y[-1][0]) / CaMKII_total]

    # print result
    print('Figure 2a K0.5 equals ' + str(get_kd(result)))

    plt.figure(1)
    plt.semilogx(data_2a[:, 0], data_2a[:, 1], 'o')
    plt.semilogx(result[:, 0], result[:, 1], '*-')


def plot_figure2b(param, data_2b):
    t0 = np.linspace(0, 60000, num=60001)
    Cai = 0.5
    ATP = 0
    CaM_total = 1E-4
    result = []

    for CaMKII_total, fraction in data_2b:
        y0 = [CaMKII_total, 0, 0, 0, 0, 0, 0]
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, param))
        result.append((CaMKII_total, (CaMKII_total - y[-1][0]) / CaM_total))

    # print result
    print('Figure 2b K0.5 equals ' + str(get_kd(result)))

    plt.figure('2b')
    plt.semilogx(data_2b[:, 0], data_2b[:, 1], 'o')
    plt.semilogx([i[0] for i in result], [i[1] for i in result], '*-')


def plot_figure2c(param, data_2c):
    t0 = np.linspace(0, 60000, num=60001)
    CaM_total = 5E-3
    CaMKII_total = 2E-4
    y0 = [CaMKII_total, 0, 0, 0, 0]
    result = []

    for Cai, fraction in data_2c:
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, param))
        result.append((Cai, (CaMKII_total - y[-1][0]) / CaMKII_total))

    # print result
    print('Figure 2c K0.5 equals ' + str(get_kd(result)))

    plt.figure(3)
    plt.semilogx(data_2c[:, 0], data_2c[:, 1], 'o')
    plt.semilogx([i[0] for i in result], [i[1] for i in result], '*-')


def plot_figure3a(param, data_3a, curve_3a):
    t0 = np.linspace(0, 15000, num=15001)
    CaMKII_total = 62E-6
    CaMKII_initial = CaMKII_total
    Cai = 0.5
    ATP = 0.1
    temperature = 0
    y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0]
    x = np.logspace(-5, -2, num=100)
    result = np.zeros((len(x), 2))

    for i, CaM_total in enumerate(x):
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param))
        result[i] = [CaM_total, 100 * (CaMKII_total - y[-1][0] - y[-1][1]) / CaMKII_total]

    result = np.array([[CaM_total, 100 * auto / result[-1][1]] for CaM_total, auto in result])

    # print result
    print('Figure 3a K0.5 equals ' + str(get_kd(result)))

    plt.figure(1)
    plt.semilogx(data_3a[:, 0], data_3a[:, 1], 'o')
    plt.semilogx(result[:, 0], result[:, 1], '-')
    plt.semilogx(x, fit_hillequation(x, *curve_3a), '--')


def plot_figure3b(param, data_3b, curve_3b):
    t0 = np.linspace(0, 6000, num=6001)
    CaMKII_total = 5E-6
    CaMKII_initial = CaMKII_total
    Cai = 0.5
    ATP = 0.25
    temperature = 30
    y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0]
    # param = np.array([2.1, 1.4E-4, 1.9E-3, 1.4E-7, 1.9E-6, 3E-4])
    x = np.logspace(-5, -2, num=100)
    result = np.zeros((len(x), 2))

    for i, CaM_total in enumerate(x):
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param))
        result[i] = [CaM_total, 100 * (CaMKII_total - y[-1][0] - y[-1][1]) / CaMKII_total]

    # print result
    print('Figure 3b K0.5 equals ' + str(get_kd(result)))

    plt.figure(2)
    plt.semilogx(data_3b[:, 0], data_3b[:, 1], 'o')
    plt.semilogx(result[:, 0], result[:, 1], '-')
    plt.semilogx(x, fit_hillequation(x, *curve_3b), '--')


def plot_figure3c(param, data_3c):
    t0 = np.linspace(0, 300000, num=300001)
    CaMKII_total = 0.2E-3
    CaMKII_initial = CaMKII_total
    CaM_total = 50E-3
    ATP = 2
    temperature = 0
    y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0]
    result = np.zeros((len(data_3c), 2))

    for i, (Cai, fraction) in enumerate(data_3c):
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param))
        result[i] = [Cai, 100 * (CaMKII_total - y[-1][0] - y[-1][1]) / CaMKII_total]

    # print result
    print('Figure 3c K0.5 equals ' + str(get_kd(result)))

    plt.figure(3)
    plt.semilogx(data_3c[:, 0], data_3c[:, 1], 'o')
    plt.semilogx(result[:, 0], result[:, 1], '*-')


def plot_figure4(param_a, param_d, data_41, curve_41, data_42, curve_42, data_44, curve_44):
    plt.figure(1, figsize=(3.24, 2.3), dpi=300)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.axes([0.15, 0.18, 0.8, 0.77])

    # plot 1Hz
    # x = np.linspace(0, 100)
    # plt.plot(x, fit_hillequation(x, *curve_41), 'b-')
    plt.plot(data_41[:, 0], data_41[:, 1], 'sb', ms=5, label='   ')
    plt.axis([0, 100, 0, 70])

    t0 = np.linspace(0, 100000, num=100001)  # simulate 100s for the 1Hz case
    CaMKII_total = 2E-5  # because the CaMKII is fixed in a tube where fluid is flowing,
    # so the concentration is unknown, this value is estimated.
    CaM_total = 1E-4
    Cai = 0.5
    ATP = 0.25
    temperature = 30
    y0 = [CaMKII_total, 0, 0, 0, 0, 0, 0]

    y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param_a, 1), hmax=10.0)
    result = [(CaMKII_total - II - II_CaM) / CaMKII_total * 100 for (II, II_CaM, IIP_CaM, d1, d2, d3, d4) in y]
    plt.plot(t0 / 1000, result, 'b--', label=r'$\mathrm{1Hz}$', lw=1.5)

    y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param_d, 1), hmax=10.0)
    result = [(CaMKII_total - II - II_CaM) / CaMKII_total * 100 for (II, II_CaM, IIP_CaM, d1, d2, d3, d4) in y]
    plt.plot(t0 / 1000, result, 'b-', label='   ', lw=1.5)

    # plot 2.5Hz
    # x = np.linspace(0, 45)
    # plt.plot(x, fit_hillequation(x, *curve_42), 'b-')
    plt.plot(data_42[:, 0], data_42[:, 1], '^g', ms=5, label='     ')

    t0 = np.linspace(0, 45000, num=45001)  # simulate 45s for the 2.5Hz case
    y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param_a, 2.5), hmax=10.0)
    result = [(CaMKII_total - II - II_CaM) / CaMKII_total * 100 for (II, II_CaM, IIP_CaM, d1, d2, d3, d4) in y]
    plt.plot(t0 / 1000, result, 'g--', label=r'$\mathrm{2.5Hz}$')

    y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param_d, 2.5), hmax=10.0)
    result = [(CaMKII_total - II - II_CaM) / CaMKII_total * 100 for (II, II_CaM, IIP_CaM, d1, d2, d3, d4) in y]
    plt.plot(t0 / 1000, result, 'g-', label='     ')

    # plot 4Hz
    plt.plot(data_44[:, 0], data_44[:, 1], 'or', ms=5, label='   ')
    # x = np.linspace(0, 30)
    # plt.plot(x, fit_hillequation(x, *curve_44), 'b-')

    t0 = np.linspace(0, 30000, num=30001)  # simulate 30s for the 4Hz case
    y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param_a, 4), hmax=10.0)
    result = [(CaMKII_total - II - II_CaM) / CaMKII_total * 100 for (II, II_CaM, IIP_CaM, d1, d2, d3, d4) in y]
    plt.plot(t0 / 1000, result, 'r--', label=r'$\mathrm{4Hz}$')

    y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param_d, 4), hmax=10.0)
    result = [(CaMKII_total - II - II_CaM) / CaMKII_total * 100 for (II, II_CaM, IIP_CaM, d1, d2, d3, d4) in y]
    plt.plot(t0 / 1000, result, 'r-', label='   ')

    # plt.legend(bbox_to_anchor=(0., 1.02, 1., 0.2), loc=3,
    #            ncol=3, mode="expand", borderaxespad=0., fontsize=9)
    plt.tick_params(axis='both', which='major', labelsize=9)
    plt.xlabel(r'$\mathrm{Time}\;(s)$', fontsize=9, labelpad=2.5)
    plt.ylabel(r'$\mathrm{Autophosphorylation}\;(\%)$', fontsize=9)
    plt.savefig('/Users/fairly/Desktop/8-15-figure4.pdf')


def plot_figure5a(param, data_5a, curve_5a):
    t0 = np.linspace(0, 60000, num=60001)
    Cai = 0.5
    CaM_total = 1E-4
    ATP = 0
    temperature = 0
    result = []

    x = np.logspace(-6, -2, num=100)
    for CaMKII_total in x:
        y0 = [CaMKII_total, 0, 0, 0, 0, 0, 0]
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param))
        result.append((CaMKII_total, (y[-1][1] + y[-1][2]) / CaM_total))

    # print result
    print('Figure 5a K0.5 equals ' + str(get_kd(result)))

    fig = plt.figure(1, figsize=(2.76, 1.9), dpi=300)
    fig.text(0.02, 0.93, 'A')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    plt.axes([0.17, 0.18, 0.74, 0.74])
    plt.semilogx(data_5a[:, 0], data_5a[:, 1] * 100, 'o', markersize=5)
    plt.semilogx([i[0] for i in result], [i[1] * 100 for i in result], '-', linewidth=1.5)
    # plt.semilogx(x, fit_hillequation(x, *curve_5a) * 100, '--', linewidth=1.5)
    plt.tick_params(axis='both', which='major', labelsize=9)
    plt.xlabel(r'$\mathrm{[CaMKII]_{total}}\;(mM)$', fontsize=9, labelpad=2)
    plt.ylabel(r'$\mathrm{CaM\;bound\;fraction}\;(\%)$', fontsize=9, labelpad=1.5)
    plt.savefig('/Users/fairly/Desktop/8-15-a.pdf')


def plot_figure5b(param, data_5b, curve_5b):
    t0 = np.linspace(0, 15000, num=15001)
    CaMKII_total = 62E-6
    CaMKII_initial = CaMKII_total
    Cai = 0.5
    ATP = 0.1
    temperature = 0
    y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0]
    x = np.logspace(-6, -2, num=100)
    result = np.zeros((len(x), 2))

    for i, CaM_total in enumerate(x):
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param))
        result[i] = [CaM_total, (CaMKII_total - y[-1][0] - y[-1][1]) / CaMKII_total]

    result = np.array([[CaM_total, auto / result[-1][1]] for CaM_total, auto in result])

    # print result
    print('Figure 5b K0.5 equals ' + str(get_kd(result)))

    fig = plt.figure(2, figsize=(2.76, 1.9), dpi=300)
    fig.text(0.02, 0.93, 'B')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    # fig.set_size_inches(3.24, 2.43)
    plt.axes([0.17, 0.18, 0.74, 0.74])
    plt.semilogx(data_5b[:, 0], data_5b[:, 1] * 100, 'o', markersize=5)
    plt.semilogx([i[0] for i in result], [i[1] * 100 for i in result], '-', linewidth=1.5)
    # plt.semilogx(x, 100 * fit_hillequation(x, *curve_5b), '--', linewidth=1.5)
    plt.tick_params(axis='both', which='major', labelsize=9)
    plt.xlabel(r'$\mathrm{[CaM]_{total}}\;(mM)$', fontsize=9, labelpad=2)
    plt.ylabel(r'$\mathrm{Autophosphorylation}\;(\%)$', fontsize=9, labelpad=1.5)
    plt.savefig('/Users/fairly/Desktop/8-15-b.pdf')


def plot_figure5c(param, data_5c, curve_5c):
    t0 = np.linspace(0, 60000, num=60001)
    CaMKII_total = 6.2E-6
    CaMKII_initial = CaMKII_total
    Cai = 0.5
    ATP = 0.1
    temperature = 30
    y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0]
    x = np.logspace(-7, -2, num=100)
    result = np.zeros((len(x), 2))

    for i, CaM_total in enumerate(x):
        y = odeint(CaMKII_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param))
        # result[i] = [CaM_total, (CaMKII_total - y[-1][0]) / CaMKII_total]
        tmp = [y[-1][0], y[-1][1], y[-1][2], CaMKII_total - y[-1][0] - y[-1][1] - y[-1][2]]
        result[i] = [CaM_total, kinase_activity(*tmp)]

    result_m = max([i for a, i in result])
    # result = np.array([[CaM_total, k_a / result[-1][1]] for CaM_total, k_a in result])
    result = np.array([[CaM_total, k_a / result_m] for CaM_total, k_a in result])

    # print result
    print('Figure 5c K0.5 equals ' + str(get_kd(result)))

    fig = plt.figure(3, figsize=(2.76, 1.9), dpi=300)
    fig.text(0.02, 0.93, 'C')
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    # fig.set_size_inches(3.24, 2.43)
    plt.axes([0.17, 0.18, 0.74, 0.74])
    plt.semilogx(data_5c[:, 0], data_5c[:, 1] * 100, 'o', markersize=5)
    plt.semilogx([i[0] for i in result], [i[1] * 100 for i in result], '-', linewidth=1.5)
    # plt.semilogx(x, fit_hillequation(x, *curve_5c) * 100, '--', linewidth=1.5)
    plt.tick_params(axis='both', which='major', labelsize=9)
    plt.ylabel(r'$\mathrm{Kinase\;Activity}\;(\%)$', fontsize=9, labelpad=1.5)
    plt.xlabel(r'$\mathrm{[CaM]_{total}}\;(mM)$', fontsize=9, labelpad=2)
    plt.savefig('/Users/fairly/Desktop/8-15-c.pdf')


def plot_ROS(param, data_ROS):
    t0 = np.linspace(0, 60000, num=60001)
    CaMKII_total = 6.2E-4
    CaMKII_initial = CaMKII_total
    CaM_total = 1E-3
    Cai = 200E-3
    ATP = 0
    temperature = 0
    y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0, 0, 0]
    x = np.logspace(-6, -3, num=100)
    MsrA = 0
    result = np.zeros((len(x), 2))

    # TODO the data used for plotting here are a little different from err_ROS
    tmp_max = max([k_a for k_a in data_ROS[:, 1]])
    data_ROS = np.array([[ROS, k_a / tmp_max * 100] for ROS, k_a in data_ROS])

    for i, ROS in enumerate(x):
        y = odeint(CaMKIIROS_transition, y0, t0, args=(
            CaM_total, CaMKII_total, Cai, ATP, temperature, ROS, MsrA, param))
        # result[i] = [ROS, (CaMKII_total - y[-1][0]) / CaMKII_total]
        # tmp = [y[-1][0], y[-1][1], y[-1][2],
        #        CaMKII_total - y[-1][0] - y[-1][1] - y[-1][2] - y[-1][3] - y[-1][4], y[-1][3], y[-1][4]]
        # result[i] = [ROS, kinase_activity(*tmp)]
        tmp = y[-1][3:5]
        result[i] = [ROS, sum(tmp)]

    result_m = max([i for a, i in result])
    # result = np.array([[ROS, k_a / result[-1][1]] for ROS, k_a in result])
    result = np.array([[ROS, k_a / result_m * 100] for ROS, k_a in result])

    # result_m = max([i for a, i in result])
    # # result = np.array([[ROS, k_a / result[-1][1]] for ROS, k_a in result])
    # result = np.array([[ROS, k_a / result_m] for ROS, k_a in result])

    # print result
    # print('Figure 5c K0.5 equals ' + str(get_kd(result)))

    fig = plt.figure(3, figsize=(2.76, 1.9), dpi=300)
    plt.rc('text', usetex=True)
    plt.rc('font', family='serif')
    # fig.set_size_inches(3.24, 2.43)
    plt.axes([0.17, 0.18, 0.74, 0.74])
    plt.semilogx(data_ROS[:, 0], data_ROS[:, 1], 'o', ms=5)
    plt.semilogx([i[0] for i in result], [i[1] for i in result], '-', linewidth=1.5)
    plt.ylabel(r'$\mathrm{Kinase\;Activity}\;(\%)$', fontsize=9, labelpad=1.5)
    plt.xlabel('$\mathrm{H_{2}O_{2}}\;(mM)$', fontsize=9, labelpad=2)
    plt.tick_params(axis='both', which='major', labelsize=9)
    plt.axis([0, 0.0010, 0, 105])
    plt.savefig('/Users/fairly/Desktop/ROS.pdf')


def set_environment(data):
    for expriment in data:
        if expriment[0] == '5a':
            t0 = np.linspace(0, 60000, num=60001)
            Cai = 0.5
            CaM_total = 1E-4
            ATP = 0
            temperature = 0
            expriment.append([t0, Cai, CaM_total, ATP, temperature])
    return data


def xls(figure_name):
    data = xlrd.open_workbook('steady-state-I-V.xlsx')
    table = data.sheet_by_name(figure_name)  # sheet_index begins from 0
    r = table.nrows  # number of rows
    c = table.ncols
    rtn = np.zeros((r, c))
    for i in range(r):
        rtn[i] = table.row_values(i)
    return rtn  # return an array just looks like the table in excel


def test_CaM():
    CaM_total = 6E-3
    y0 = [0, 0, 0, 0]
    t0 = np.linspace(0, 10000, num=1000001)

    y = odeint(CaM_transition, y0, t0, args=(CaM_total,))

    fig = plt.figure()

    ax1 = fig.add_subplot(2, 1, 1)
    ax1.plot(t0, [get_Cai(t) for t in t0])

    ax2 = fig.add_subplot(2, 1, 2)
    ax2.plot(t0, y[:, 1], label='CaMCa2')
    ax2.set_ylim(0, 6E-4)
    ax2.legend()

    ax3 = ax2.twinx()
    ax3.plot(t0, y[:, 3], label='CaMCa4')
    ax3.set_ylim(0, 3E-6)
    ax3.legend()

    plt.show()


if __name__ == '__main__':
    # print_k()
    # test_CaM()
    # data = xls('2b')
    # fit_figure4()
    fit_figure5()
    # fit_ROS()
    # curvfit(data)
    # investigation()
    # y0 = [CaMKII_initial, 0, 0, 0, 0, 0, 0]
    # y = odeint(CaMKIIROS_transition, y0, t0, args=(CaM_total, CaMKII_total, Cai, ATP, temperature, param))
    pass
