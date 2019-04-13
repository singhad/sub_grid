import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D

def make_pdf(s, s_bar, sigma_s):
    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_bar)/sigma_s)**2)))
    return pdf

def calc_integral(s, pdf, n_H_mean, X_H2, ds):
    integ = 0.0
    for i in range(1000):
        integ += (np.exp(s[i])*pdf[i]*X_H2[i]*ds)
    return integ

def calc_lambda_jeans(n_H, c_s, G, K_b, m_p):
    lambda_jeans = c_s / np.sqrt(4* np.pi * G * n_H * m_p)
    return lambda_jeans

def calc_n_LW(n_H, G_o, lambda_jeans, m_p):
    kappa = 1000 * m_p
    rad_field_outside = G_o #in solar units
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    n_LW = rad_field_outside * exp_tau
    return n_LW

def calc_n_LW_ss(n_H, n_H2, G_o, lambda_jeans, m_p):
    kappa = 1000 * m_p
    rad_field_outside = G_o #in solar units
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    N_H2 = n_H2*lambda_jeans
    term1 = (0.965/((1+(N_H2/5e14))**2))
    term2 = ( (0.035/np.sqrt(1+(N_H2/5e14))) * np.exp(-1*np.sqrt(1+(N_H2/5e14))/1180) )
    S_H2 = term1 + term2
    n_LW_ss = rad_field_outside * exp_tau * S_H2
    return n_LW_ss, S_H2, N_H2

def calc_X_H2(n_H, Z, n_LW):
    DC = 1.7e-11
    CC = 2.5e-17            #cm3 s-1
    numerator = DC * n_LW
    denominator = CC * Z * n_H
    X_H2 = 1 / (2 + (numerator/denominator) )
    return X_H2

def self_shielding_iterations(n_H, G_o, lambda_jeans, Z, m_p):
    X_H2 = np.zeros(1000)
    n_LW = np.zeros(1000)
    n_H2 = np.zeros(1000)
    n_LW_ss = np.zeros(1000)
    S_H2_ss = np.zeros(1000)
    N_H2_ss = np.zeros(1000)
    X_H2_ss = np.zeros(1000)
    n_H2_ss = np.zeros(1000)
    ctr = 15
    i = 0
    n_LW = calc_n_LW(n_H, G_o, lambda_jeans, m_p)
    X_H2 = calc_X_H2(n_H, Z, n_LW)
    n_H2 = n_H * X_H2
    n_H2_ss = n_H2
    while i<ctr:
        n_LW_ss, S_H2_ss, N_H2_ss = calc_n_LW_ss(n_H, n_H2_ss, G_o, lambda_jeans, m_p)
        X_H2_ss = calc_X_H2(n_H, Z, n_LW_ss)
        n_H2_ss = n_H * X_H2_ss
        i += 1
    return n_LW, n_LW_ss, S_H2_ss, N_H2_ss, X_H2_ss, n_H2_ss, n_H2

def calc_X_CO(n_H, n_H2, n_LW):
    rate_CHX = 5.e-10 * n_LW
    rate_CO = 1.e-10 * n_LW
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

def inside_loop(M, n_H_mean, Z, G_o):
    T = 10 #K
    m_p = 1.672621777e-24   # g
    T_mean = 10.            #K
    K_b = 1.38064852e-16    # ergs K-1
    G = 6.67408e-8          # dyne cm^2 g^-2
    c_s = np.sqrt(K_b * T / m_p)
    s = np.zeros(1000)
    pdf = np.zeros(1000)
    lambda_jeans = np.zeros(1000)
    n_H = np.zeros(1000)
    X_CO = np.zeros(1000)
    n_CO = np.zeros(1000)
    integral1 = 0
    sigma_s = np.sqrt(np.log(1 + ((0.3 * M)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -7*sigma_s + s_bar
    smax = 7*sigma_s + s_bar
    ds = (smax - smin)/1000
    for i in range(1000):
        s[i] = smin + i*ds
    pdf = make_pdf(s, s_bar, sigma_s)
    n_H = n_H_mean * np.exp(s)
    lambda_jeans = calc_lambda_jeans(n_H, c_s, G, K_b, m_p)
    n_LW, n_LW_ss, S_H2_ss, N_H2_ss, X_H2_ss, n_H2_ss, n_H2 = self_shielding_iterations(n_H, G_o, lambda_jeans, Z, m_p)
    integral1 = calc_integral(s, pdf, n_H_mean, X_H2_ss, ds)
    X_H2_bar = integral1
    X_CO = calc_X_CO(n_H, n_H2, n_LW)
    n_CO = calc_n_CO(n_H, X_CO)
    return X_H2_bar

def varying_M():
    fig, ax = plt.subplots()
    Z = 1
    G_o = 1
    mach_no_arr = np.logspace(-2, 2, 5)
    n_H_mean_arr = np.logspace(-2, 3.5, 40)
    X_H2_bar = np.zeros(len(n_H_mean_arr))
    label = "M "
    color_arr = ['r','g','b','y','k']
    custom_lines = [Line2D([0], [0], color='r', lw=4),
                    Line2D([0], [0], color='g', lw=4),
                    Line2D([0], [0], color='b', lw=4),
                    Line2D([0], [0], color='y', lw=4),
                    Line2D([0], [0], color='k', lw=4)]
    for m in range(len(mach_no_arr)):
        mach_no = mach_no_arr[m]
        color = str(color_arr[m])
        for y in range(len(n_H_mean_arr)):
            n_H_mean = n_H_mean_arr[y]
            X_H2_bar = inside_loop(mach_no, n_H_mean, Z, G_o)
            plt.scatter(np.log10(n_H_mean), X_H2_bar, color=color)

    ax.legend(  custom_lines,
                [   label + '= 1e-2',
                    label + '= 1e-1',
                    label + '= 1e0',
                    label + '= 1e+1',
                    label + '= 1e+2'  ],
                loc = 'lower right'
                    )
    plt.xlabel('log(n_H_mean)')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H_mean) vs X_H2_bar - M=varied, Z=1, G_o=1')
    plt.savefig(os.path.join('log(n_H_mean)vsX_H2_bar--M.png'.format()))
    plt.clf()
    return X_H2_bar

def varying_Z():
    fig, ax = plt.subplots()
    mach_no = 10
    G_o = 1
    Z_arr = np.array([0.001, 0.01, 0.1, 1.])
    n_H_mean_arr = np.logspace(-2, 3.5, 40)
    X_H2_bar = np.zeros(len(n_H_mean_arr))
    label = "Z "
    color_arr = ['r','g','b','y']
    custom_lines = [Line2D([0], [0], color='r', lw=4),
                    Line2D([0], [0], color='g', lw=4),
                    Line2D([0], [0], color='b', lw=4),
                    Line2D([0], [0], color='y', lw=4)]
    for z in range(0, len(Z_arr)):
        Z = Z_arr[z]
        color = str(color_arr[z])
        for y in range(len(n_H_mean_arr)):
            n_H_mean = n_H_mean_arr[y]
            X_H2_bar = inside_loop(mach_no, n_H_mean, Z, G_o)
            plt.scatter(np.log10(n_H_mean), X_H2_bar, color=color)
    plt.xlabel('log(n_H_mean)')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H_mean) vs X_H2_bar - M=10, Z=varied, G_o=1')
    ax.legend(  custom_lines,
                [   label + '= 0.001',
                    label + '= 0.010',
                    label + '= 0.100',
                    label + '= 1.000'  ],
                loc = 'lower right'
                    )
    plt.savefig(os.path.join('log(n_H_mean)vsX_H2_bar--Z.png'.format()))
    plt.clf()
    return X_H2_bar

def varying_G_o():
    fig, ax = plt.subplots()
    mach_no = 10
    Z = 1
    G_o_arr = np.array([1., 10., 50., 100.])
    n_H_mean_arr = np.logspace(-2, 3.5, 40)
    X_H2_bar = np.zeros(len(n_H_mean_arr))
    label = "G_o "
    color_arr = ['r','g','b','y']
    custom_lines = [Line2D([0], [0], color='r', lw=4),
                    Line2D([0], [0], color='g', lw=4),
                    Line2D([0], [0], color='b', lw=4),
                    Line2D([0], [0], color='y', lw=4)]
    for g in range(0, len(G_o_arr)):
        G_o = G_o_arr[g]
        color = str(color_arr[g])
        for y in range(len(n_H_mean_arr)):
            n_H_mean = n_H_mean_arr[y]
            X_H2_bar = inside_loop(mach_no, n_H_mean, Z, G_o)
            plt.scatter(np.log10(n_H_mean), X_H2_bar, color=color)
    plt.xlabel('log(n_H_mean)')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H_mean) vs X_H2_bar - M=10, Z=1, G_o=varied')
    ax.legend(  custom_lines,
                [   label + '= 1',
                    label + '= 10',
                    label + '= 50',
                    label + '= 100'  ],
                loc = 'lower right'
                    )
    plt.savefig(os.path.join('log(n_H_mean)vsX_H2_bar--G_o.png'.format()))
    plt.clf()
    return X_H2_bar

if __name__=='__main__':
    path = 'for X_H2_bar_self_shielding'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)

    # order of variables:
    # s, smin, smax, sigma_s, n_H, lambda_jeans, X_H2, pdf, integral1, X_H2_bar
    #g1, g2, g3, g4, g5, g6, g7, g8, g9,
    g10 = varying_G_o()  #varying G_o
    #m1, m2, m3, m4, m5, m6, m7, m8, m9,
    m10 = varying_M()  #varying mach_no
    #z1, z2, z3, z4, z5, z6, z7, z8, z9,
    z10 = varying_Z()  #varying metallicity

'''
n_LW[i] = calc_n_LW(n_H[i], G_o, lambda_jeans[i])
    X_H2_a[i] = calc_X_H2(n_H[i], Z, n_LW[i])
    n_H2_a[i] = n_H[i] * X_H2_a[i]

    n_LW_1[i], S_H2_1[i], N_H2_1[i] = calc_n_LW_ss(n_H[i], n_H2_a[i], G_o, lambda_jeans[i])
    X_H2_1[i] = calc_X_H2(n_H[i], Z, n_LW_1[i])
    n_H2_1[i] = n_H[i] * X_H2_1[i]

    n_LW_2[i], S_H2_2[i], N_H2_2[i] = calc_n_LW_ss(n_H[i], n_H2_1[i], G_o, lambda_jeans[i])
    X_H2_2[i] = calc_X_H2(n_H[i], Z, n_LW_2[i])
    n_H2_2[i] = n_H[i] * X_H2_2[i]

    n_LW_3[i], S_H2_3[i], N_H2_3[i] = calc_n_LW_ss(n_H[i], n_H2_2[i], G_o, lambda_jeans[i])
    X_H2_3[i] = calc_X_H2(n_H[i], Z, n_LW_3[i])
    n_H2_3[i] = n_H[i] * X_H2_3[i]

    n_LW_4[i], S_H2_4[i], N_H2_4[i] = calc_n_LW_ss(n_H[i], n_H2_3[i], G_o, lambda_jeans[i])
    X_H2_4[i] = calc_X_H2(n_H[i], Z, n_LW_4[i])
    n_H2_4[i] = n_H[i] * X_H2_4[i]

    n_LW_5[i], S_H2_5[i], N_H2_5[i] = calc_n_LW_ss(n_H[i], n_H2_4[i], G_o, lambda_jeans[i])
    X_H2_5[i] = calc_X_H2(n_H[i], Z, n_LW_5[i])
    n_H2_5[i] = n_H[i] * X_H2_5[i]

    n_LW_6[i], S_H2_6[i], N_H2_6[i] = calc_n_LW_ss(n_H[i], n_H2_5[i], G_o, lambda_jeans[i])
    X_H2_6[i] = calc_X_H2(n_H[i], Z, n_LW_6[i])
    n_H2_6[i] = n_H[i] * X_H2_6[i]

    n_LW_7[i], S_H2_7[i], N_H2_7[i] = calc_n_LW_ss(n_H[i], n_H2_6[i], G_o, lambda_jeans[i])
    X_H2_7[i] = calc_X_H2(n_H[i], Z, n_LW_7[i])
    n_H2_7[i] = n_H[i] * X_H2_7[i]

    n_LW_8[i], S_H2_8[i], N_H2_8[i] = calc_n_LW_ss(n_H[i], n_H2_7[i], G_o, lambda_jeans[i])
    X_H2_8[i] = calc_X_H2(n_H[i], Z, n_LW_8[i])
    n_H2_8[i] = n_H[i] * X_H2_8[i]

    n_LW_9[i], S_H2_9[i], N_H2_9[i] = calc_n_LW_ss(n_H[i], n_H2_8[i], G_o, lambda_jeans[i])
    X_H2_9[i] = calc_X_H2(n_H[i], Z, n_LW_9[i])
    n_H2_9[i] = n_H[i] * X_H2_9[i]

    n_LW_10[i], S_H2_10[i], N_H2_10[i] = calc_n_LW_ss(n_H[i], n_H2_9[i], G_o, lambda_jeans[i])
    X_H2_10[i] = calc_X_H2(n_H[i], Z, n_LW_10[i])
    n_H2_10[i] = n_H[i] * X_H2_10[i]

    n_LW_ss[i], S_H2[i], N_H2[i] = calc_n_LW_ss(n_H[i], n_H2_10[i], G_o, lambda_jeans[i])
    X_H2[i] = calc_X_H2(n_H[i], Z, n_LW_ss[i])
    n_H2[i] = n_H[i] * X_H2[i]
'''
