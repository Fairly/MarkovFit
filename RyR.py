__author__ = 'fairly'

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import odeint


def CaMKIIRYR_transition(y, t0, CaM_total, CaMKII_total, Cai, ATP, temperature, ROS, MsrA,
                         param):
    CaMKII, CaMKII_CaMCa4, CaMKIIP_CaMCa4, OX, OXA, CaMCa, CaMCa2, CaMCa3, CaMCa4, RyRr, RyRo, RyRi, Caryr, Csqn_b = y
    CaMKIIP = CaMKII_total - CaMKII - CaMKII_CaMCa4 - CaMKIIP_CaMCa4 - OX - OXA
    CaM = CaM_total - CaMCa - CaMCa2 - CaMCa3 - CaMCa4 - CaMKII_CaMCa4 - CaMKIIP_CaMCa4 - OX

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

    kcat, k_asso, k_disso, k_dissoCa, k_ox = [2.73e-05, 1.15e+00, 1e-08, 2e-03, 1.47e-05]
    # kcat = args
    # k_asso = 2.1
    # k_disso = 0.7E-4
    # k_dissoCa = 0.95E-3
    k_disso2 = k_disso / 1000
    k_dissoCa2 = k_dissoCa / 1000
    # k_ox = 1.952869E-6
    k_red = 0.01432

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
    d4 = -CA4 - A1 + k_disso * cons_cam * CaMKII_CaMCa4 - C2 + k_disso2 * cons_cam * CaMKIIP_CaMCa4 + k_disso * \
                                                                                                      cons_cam * OX - F2

    cellRadius = 10.25
    cellLength = 100.0
    Vcell = 3.14159265358979 * pow(cellRadius, 2.0) * cellLength * 1.0e-15
    Vmyo = 0.65 * Vcell
    Vsr = 0.035 * Vcell     # geometry

    kon_csqn = 100.0
    Bmax_Csqn = 140.0e-3 * Vmyo / Vsr
    koff_csqn = 65.0
    MaxSR = 15.0
    MinSR = 1.0
    ec50SR = 0.45
    kom = 0.06
    kim = 0.005
    ks = 25.0
    koCa = 10.0
    kiCa = 0.5
    kCaSR = MaxSR - (MaxSR - MinSR) / (1.0 + pow(ec50SR / Caryr, 2.5))
    koSRCa = koCa / kCaSR
    kiSRCa = kiCa * kCaSR
    ec50SR, kom, kim, ks, koCa, kiCa = param    # redefine some parameters, an easy-to-use mothed

    RI = 1.0 - RyRr - RyRo - RyRi
    dRyRr = kim * RI - kiSRCa * Cai * RyRr \
            - (koSRCa * pow(Cai, 2.0) * RyRr - kom * RyRo)
    dRyRo = koSRCa * pow(Cai, 2.0) * RyRr - kom * RyRo \
            - (kiSRCa * Cai * RyRo - kim * RyRi)
    dRyRi = kiSRCa * Cai * RyRo - kim * RyRi \
            - (kom * RyRi - koSRCa * pow(Cai, 2.0) * RI)
    J_SRCarel = ks * RyRo / 1.0 * (Caryr - Cai)
    dCaryr = J_SRCarel - (kon_csqn * Caryr * (Bmax_Csqn - Csqn_b) - koff_csqn * Csqn_b)
    dCsqn_b = kon_csqn * dCaryr * (Bmax_Csqn - Csqn_b) - koff_csqn * Csqn_b

    return [dCaMKII, dCaMKII_CaMCa4, dCaMKIIP_CaMCa4, dOX, dOXA, d1, d2, d3, d4, dRyRr, dRyRo, dRyRi, dCaryr, dCsqn_b]


if __name__ == '__main__':
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

    param = [0.45, 0.06, 0.005, 25, 10, 0.5]

    tmp_max = max([k_a for k_a in data_ROS[:, 1]])
    data_ROS = np.array([[ROS, k_a / tmp_max * 100] for ROS, k_a in data_ROS])

    for i, ROS in enumerate(x):
        y = odeint(CaMKIIRYR_transition, y0, t0, args=(
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

    plt.figure(1)
    plt.plot(data_ROS[:, 0], data_ROS[:, 1], 'o')
    plt.plot([i[0] for i in result], [i[1] for i in result], '-')
    plt.ylabel('Kinase Activity (%)')
    plt.xlabel('H$_{2}$O$_{2}$ (mM)')
    plt.axis([0, 0.0010, 0, 105])
