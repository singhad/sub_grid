import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D

def make_pdf(s, s_bar, sigma_s):
    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_bar)/sigma_s)**2)))
    return pdf

def calc_lambda_jeans(n_H):
    m_p = 1.672621777e-24   # g
    T_mean = 10.            #K
    K_b = 1.38064852e-16    # ergs K-1
    G = 6.67408e-8          # dyne cm^2 g^-2
    lambda_jeans = ((np.sqrt(K_b * T_mean / m_p)) / np.sqrt(4* np.pi * G * n_H * m_p))
    return lambda_jeans

def calc_n_LW(n_H, G_o, lambda_jeans):
    m_p = 1.672621777e-24   # g
    kappa = 1000 * m_p
    rad_field_outside = G_o #in solar units
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    n_LW = rad_field_outside * exp_tau
    return n_LW

def calc_n_LW2(n_H, n_H2, G_o, lambda_jeans):
    m_p = 1.672621777e-24   # g
    kappa = 1000 * m_p
    rad_field_outside = G_o #in solar units
    exp_tau_H = np.exp(-kappa * n_H * lambda_jeans)
    exp_tau_H2 = np.exp(-kappa * n_H2 * lambda_jeans)
    n_LW2 = rad_field_outside * exp_tau_H * exp_tau_H2
    return n_LW2

def calc_X_H2(n_H, Z, n_LW):
    DC = 1.7e-11
    CC = 2.5e-17            #cm3 s-1
    numerator = DC * n_LW
    denominator = CC * Z * n_H
    X_H2 = 1 / (2 + (numerator/denominator) )
    return X_H2

def calc_X_CO(n_H, n_H2, n_LW2):
    rate_CHX = 5.e-10 * n_LW2
    rate_CO = 1.e-10 * n_LW2
    x0 = 2.e-4
    k0 = 5.e-16 #cm3 s-1
    k1 = 5.e-10 #cm3 s-1
    factor_beta = rate_CHX/(n_H*k1*x0)
    beta = 1./(1.+factor_beta)
    factor_CO = rate_CO/(n_H2*k0*beta)
    X_CO = 1./(1.+factor_CO)
    return X_CO

def calc_n_CO(n_H, X_CO):
    abundance_Ctot = 1e-4 # n_C/n_H as defined by nucleosynthesis
    return n_H * abundance_Ctot * X_CO # CO/cc

def make_integral2():
    s = np.zeros(1000)
    pdf = np.zeros(1000)
    integral2 = 0
    mach_no = 5
    sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -4*sigma_s + s_bar
    smax = 4*sigma_s + s_bar
    ds = (smax - smin)/1000
    n_H_mean = 1e4
    Z = 1
    G_o = 1
    for i in range(0, 1000):
        s[i] = smin + i*ds
        pdf = make_pdf(s, s_bar, sigma_s)
        integral2 += np.exp(s[i]) * pdf[i] * ds   #this should be ~1
    #plotting(n_H, pdf, lambda_jeans, X_H2)
    return integral2


def varying_M():
    fig, ax = plt.subplots()
    Z = 1
    G_o = 1
    mach_no_arr = np.array([1., 5., 10., 50.])
    n_H_mean_arr = np.logspace(1, 5, 40)
    label = "M "
    color_arr = ['r','g','b','y']
    custom_lines = [Line2D([0], [0], color='r', lw=4),
                    Line2D([0], [0], color='g', lw=4),
                    Line2D([0], [0], color='b', lw=4),
                    Line2D([0], [0], color='y', lw=4)]
    for m in range(0, len(mach_no_arr)):
        mach_no = mach_no_arr[m]
        color = str(color_arr[m])
        sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
        s_bar = -0.5*(sigma_s**2)
        smin = -4*sigma_s + s_bar
        smax = 4*sigma_s + s_bar
        ds = (smax - smin)/1000
        value = mach_no
        for y in range(0, len(n_H_mean_arr)):
            n_H_mean = n_H_mean_arr[y]
            s = np.zeros(1000)
            pdf = np.zeros(1000)
            lambda_jeans = np.zeros(1000)
            X_H2 = np.zeros(1000)
            n_H = np.zeros(1000)
            n_LW = np.zeros(1000)
            n_LW2 = np.zeros(1000)
            n_H2 = np.zeros(1000)
            X_CO = np.zeros(1000)
            n_CO = np.zeros(1000)
            integral2 = 0
            X_CO_bar = np.zeros(1000)
            for i in range(0, 1000):
                s[i] = smin + i*ds
                n_H[i] = n_H_mean * np.exp(s[i])
                lambda_jeans[i] = calc_lambda_jeans(n_H[i])
                n_LW[i] = calc_n_LW(n_H[i], G_o, lambda_jeans[i])
                X_H2[i] = calc_X_H2(n_H[i], Z, n_LW[i])
                pdf[i] = make_pdf(s[i], s_bar, sigma_s)
                n_H2[i] = n_H[i] * X_H2[i]
                n_LW2[i] = calc_n_LW2(n_H[i], n_H2[i], G_o, lambda_jeans[i])
                X_CO[i] = calc_X_CO(n_H[i], n_H2[i], n_LW2[i])
                n_CO[i] = calc_n_CO(n_H[i], X_CO[i])
            integral2 = make_integral2()
            X_CO_bar = integral2 * X_CO
            plt.scatter(np.log10(n_H_mean), X_CO_bar[i], color=color)
    plt.xlabel('log(n_H_mean)')
    plt.ylabel('X_CO_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H_mean) vs X_CO_bar - M=varied, Z=1, G_o=1')
    ax.legend(  custom_lines,
                [   label + '= 1',
                    label + '= 5',
                    label + '= 10',
                    label + '= 50'  ],
                loc = 'lower right'
                    )
    plt.savefig(os.path.join('log(n_H_mean)vsX_CO_bar--M.png'.format()))
    plt.clf()
    return s, smin, smax, sigma_s, n_H, lambda_jeans, X_H2, pdf, n_H2, integral2, X_CO, n_CO, X_CO_bar

def varying_Z():
    fig, ax = plt.subplots()
    mach_no = 5
    G_o = 1
    Z_arr = np.array([0.001, 0.01, 0.1, 1.])
    n_H_mean_arr = np.logspace(1, 5, 40)
    label = "Z "
    color_arr = ['r','g','b','y']
    custom_lines = [Line2D([0], [0], color='r', lw=4),
                    Line2D([0], [0], color='g', lw=4),
                    Line2D([0], [0], color='b', lw=4),
                    Line2D([0], [0], color='y', lw=4)]
    for z in range(0, len(Z_arr)):
        Z = Z_arr[z]
        color = str(color_arr[z])
        sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
        s_bar = -0.5*(sigma_s**2)
        smin = -4*sigma_s + s_bar
        smax = 4*sigma_s + s_bar
        ds = (smax - smin)/1000
        value = Z
        for y in range(0, len(n_H_mean_arr)):
            n_H_mean = n_H_mean_arr[y]
            s = np.zeros(1000)
            pdf = np.zeros(1000)
            lambda_jeans = np.zeros(1000)
            X_H2 = np.zeros(1000)
            n_H = np.zeros(1000)
            n_LW = np.zeros(1000)
            n_LW2 = np.zeros(1000)
            n_H2 = np.zeros(1000)
            X_CO = np.zeros(1000)
            n_CO = np.zeros(1000)
            integral2 = 0
            X_CO_bar = np.zeros(1000)
            for i in range(0, 1000):
                s[i] = smin + i*ds
                n_H[i] = n_H_mean * np.exp(s[i])
                lambda_jeans[i] = calc_lambda_jeans(n_H[i])
                n_LW[i] = calc_n_LW(n_H[i], G_o, lambda_jeans[i])
                X_H2[i] = calc_X_H2(n_H[i], Z, n_LW[i])
                pdf[i] = make_pdf(s[i], s_bar, sigma_s)
                n_H2[i] = n_H[i] * X_H2[i]
                n_LW2[i] = calc_n_LW2(n_H[i], n_H2[i], G_o, lambda_jeans[i])
                X_CO[i] = calc_X_CO(n_H[i], n_H2[i], n_LW2[i])
                n_CO[i] = calc_n_CO(n_H[i], X_CO[i])
            integral2 = make_integral2()
            X_CO_bar = integral2 * X_CO
            plt.scatter(np.log10(n_H_mean), X_CO_bar[i], color=color)
    plt.xlabel('log(n_H_mean)')
    plt.ylabel('X_CO_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H_mean) vs X_CO_bar - M=5, Z=varied, G_o=1')
    ax.legend(  custom_lines,
                [   label + '= 0.001',
                    label + '= 0.010',
                    label + '= 0.100',
                    label + '= 1.000'  ],
                loc = 'lower right'
                    )
    plt.savefig(os.path.join('log(n_H_mean)vsX_CO_bar--Z.png'.format()))
    plt.clf()
    return s, smin, smax, sigma_s, n_H, lambda_jeans, X_H2, pdf, n_H2, integral2, X_CO, n_CO, X_CO_bar

def varying_G_o():
    fig, ax = plt.subplots()
    mach_no = 5
    Z = 1
    G_o_arr = np.array([1., 10., 50., 100.])
    n_H_mean_arr = np.logspace(1, 5, 40)
    label = "G_o "
    color_arr = ['r','g','b','y']
    custom_lines = [Line2D([0], [0], color='r', lw=4),
                    Line2D([0], [0], color='g', lw=4),
                    Line2D([0], [0], color='b', lw=4),
                    Line2D([0], [0], color='y', lw=4)]
    for g in range(0, len(G_o_arr)):
        G_o = G_o_arr[g]
        color = str(color_arr[g])
        sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
        s_bar = -0.5*(sigma_s**2)
        smin = -4*sigma_s + s_bar
        smax = 4*sigma_s + s_bar
        ds = (smax - smin)/1000
        value = G_o
        for y in range(0, len(n_H_mean_arr)):
            n_H_mean = n_H_mean_arr[y]
            s = np.zeros(1000)
            pdf = np.zeros(1000)
            lambda_jeans = np.zeros(1000)
            X_H2 = np.zeros(1000)
            n_H = np.zeros(1000)
            n_LW = np.zeros(1000)
            n_LW2 = np.zeros(1000)
            n_H2 = np.zeros(1000)
            X_CO = np.zeros(1000)
            n_CO = np.zeros(1000)
            integral2 = 0
            X_CO_bar = np.zeros(1000)
            for i in range(0, 1000):
                s[i] = smin + i*ds
                n_H[i] = n_H_mean * np.exp(s[i])
                lambda_jeans[i] = calc_lambda_jeans(n_H[i])
                n_LW[i] = calc_n_LW(n_H[i], G_o, lambda_jeans[i])
                X_H2[i] = calc_X_H2(n_H[i], Z, n_LW[i])
                pdf[i] = make_pdf(s[i], s_bar, sigma_s)
                n_H2[i] = n_H[i] * X_H2[i]
                n_LW2[i] = calc_n_LW2(n_H[i], n_H2[i], G_o, lambda_jeans[i])
                X_CO[i] = calc_X_CO(n_H[i], n_H2[i], n_LW2[i])
                n_CO[i] = calc_n_CO(n_H[i], X_CO[i])
            integral2 = make_integral2()
            X_CO_bar = integral2 * X_CO
            plt.scatter(np.log10(n_H_mean), X_CO_bar[i], color=color)
    plt.xlabel('log(n_H_mean)')
    plt.ylabel('X_CO_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H_mean) vs X_CO_bar - M=5, Z=1, G_o=varied')
    ax.legend(  custom_lines,
                [   label + '= 1',
                    label + '= 10',
                    label + '= 50',
                    label + '= 100'  ],
                loc = 'lower right'
                    )
    plt.savefig(os.path.join('log(n_H_mean)vsX_CO_bar--G_o.png'.format()))
    plt.clf()
    return s, smin, smax, sigma_s, n_H, lambda_jeans, X_H2, pdf, n_H2, integral2, X_CO, n_CO, X_CO_bar

if __name__=='__main__':
    path = 'for X_CO_bar_self_shielding'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)

    # order of variables:
    # s, smin, smax, sigma_s, n_H, lambda_jeans, X_H2, pdf, n_H2, integral2, X_CO, n_CO, X_CO_bar
    g1, g2, g3, g4, g5, g6, g7, g8, g9, g10, g11, g12, g13 = varying_G_o()  #varying G_o
    m1, m2, m3, m4, m5, m6, m7, m8, m9, m10, m11, m12, m13 = varying_M()  #varying mach_no
    z1, z2, z3, z4, z5, z6, z7, z8, z9, z10, z11, z12, z13 = varying_Z()  #varying metallicity
