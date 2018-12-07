import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D

def make_pdf(s, s_bar, sigma_s):
    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_bar)/sigma_s)**2)))
    return pdf

def make_lambda_jeans(n_H, s):
    m_p = 1.672621777e-24   # g
    T_mean = 10.            #K
    K_b = 1.38064852e-16    # ergs K-1
    G = 6.67408e-8          # dyne cm^2 g^-2
    lambda_jeans = ((np.sqrt(K_b * T_mean / m_p)) / np.sqrt(4* np.pi * G * n_H * m_p))
    return lambda_jeans

def make_X_H2(n_H, Z, G_o):
    m_p = 1.672621777e-24   # g
    T_mean = 10.            #K
    K_b = 1.38064852e-16    # ergs K-1
    G = 6.67408e-8          # dyne cm^2 g^-2
    kappa = 1000 * m_p
    DC = 1.7e-11
    CC = 2.5e-17            #cm3 s-1
    rad_field_outside = G_o #in solar units
    exp_tau = np.exp(-kappa * n_H * ((np.sqrt(K_b * T_mean / m_p)) / np.sqrt(4* np.pi * G * n_H * m_p)))
    numerator = DC * rad_field_outside * exp_tau
    denominator = CC * Z * n_H
    X_H2 = 1 / (2 + (numerator/denominator) )
    return X_H2

def make_integral1():
    s = np.zeros(1000)
    pdf = np.zeros(1000)
    lambda_jeans = np.zeros(1000)
    X_H2 = np.zeros(1000)
    n_H = np.zeros(1000)
    integral1 = 0
    mach_no = 5
    sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -3*sigma_s + s_bar
    smax = 3*sigma_s + s_bar
    ds = (smax - smin)/1000
    n_H_mean = 1e5
    Z = 1
    G_o = 1
    for i in range(0, 1000):
        s[i] = smin + i*ds
        n_H = n_H_mean * np.exp(s)
        lambda_jeans = make_lambda_jeans(n_H, s)
        X_H2 = make_X_H2(n_H, Z, G_o)
        pdf = make_pdf(s, s_bar, sigma_s)
        integral1 += np.exp(s[i]) * pdf[i] * ds   #this should be ~1
    #plotting(n_H, pdf, lambda_jeans, X_H2)
    return integral1

def plotting(n_H, pdf, lambda_jeans, X_H2):
    plt.plot(np.log10(n_H), pdf)
    plt.xlabel('log(n_H)')
    plt.ylabel('pdf')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs pdf')
    plt.savefig(os.path.join('log(n_H)vspdf.png'.format()))
    plt.clf()

    plt.plot(np.log10(n_H), np.log10(pdf))
    plt.xlabel('log(n_H)')
    plt.ylabel('log(pdf)')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs log(pdf)')
    plt.savefig(os.path.join('log(n_H)vslog(pdf).png'.format()))
    plt.clf()

    plt.plot(np.log10(lambda_jeans), np.log10(n_H))
    plt.xlabel('log(lambda_jeans)')
    plt.ylabel('log(n_H)')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(lambda_jeans) vs log(n_H)')
    plt.savefig(os.path.join('log(lambda_jeans)vslog(n_H).png'.format()))
    plt.clf()

    plt.plot(np.log10(n_H), X_H2)
    plt.xlabel('log(n_H)')
    plt.ylabel('X_H2')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs X_H2')
    plt.savefig('log(n_H)vsX_H2.png'.format())
    plt.clf()

def varying_M():
    fig, ax = plt.subplots()
    Z = 1
    G_o = 1
    mach_no_arr = np.array([1., 10., 100.])
    n_H_mean_arr = np.array([1e1, 1e2, 2.5e2, 5e2, 7.5e2, 1e3, 2.5e3, 5e3, 1e4, 1e5, 1e6, 1e7])
    label = "M "
    color_arr = ['r','g','b']
    custom_lines = [Line2D([0], [0], color='r', lw=4),
                    Line2D([0], [0], color='g', lw=4),
                    Line2D([0], [0], color='b', lw=4)]
    for m in range(0, len(mach_no_arr)):
        mach_no = mach_no_arr[m]
        color = str(color_arr[m])
        sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
        s_bar = -0.5*(sigma_s**2)
        smin = -3*sigma_s + s_bar
        smax = 3*sigma_s + s_bar
        ds = (smax - smin)/1000
        value = mach_no
        for y in range(0, len(n_H_mean_arr)):
            n_H_mean = n_H_mean_arr[y]
            s = np.zeros(1000)
            pdf = np.zeros(1000)
            lambda_jeans = np.zeros(1000)
            X_H2 = np.zeros(1000)
            n_H = np.zeros(1000)
            integral1 = 0
            X_H2_bar = np.zeros(1000)
            for i in range(0, 1000):
                s[i] = smin + i*ds
                n_H[i] = n_H_mean * np.exp(s[i])
                lambda_jeans[i] = make_lambda_jeans(n_H[i], s[i])
                X_H2[i] = make_X_H2(n_H[i], Z, G_o)
                pdf[i] = make_pdf(s[i], s_bar, sigma_s)
            integral1 = make_integral1()
            X_H2_bar = integral1 * X_H2
            plt.scatter(np.log10(n_H_mean), X_H2_bar[i], color=color)
    plt.xlabel('log(n_H_mean)')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H_mean) vs X_H2_bar - varying M')
    ax.legend(  custom_lines,
                [   label + '= 1',
                    label + '= 10',
                    label + '= 100'  ],
                loc = 'lower right'
                    )
    plt.savefig(os.path.join('log(n_H_mean)vsX_H2_bar--M.png'.format()))
    plt.clf()
    return s, smin, smax, sigma_s, n_H, lambda_jeans, X_H2, pdf, integral1, X_H2_bar

def varying_Z():
    fig, ax = plt.subplots()
    mach_no = 5
    G_o = 1
    Z_arr = np.array([0.001, 0.01, 1.])
    n_H_mean_arr = np.array([1e1, 1e2, 2.5e2, 5e2, 7.5e2, 1e3, 2.5e3, 5e3, 1e4, 1e5, 1e6, 1e7])
    label = "Z "
    color_arr = ['r','g','b']
    custom_lines = [Line2D([0], [0], color='r', lw=4),
                    Line2D([0], [0], color='g', lw=4),
                    Line2D([0], [0], color='b', lw=4)]
    for z in range(0, len(Z_arr)):
        Z = Z_arr[z]
        color = str(color_arr[z])
        sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
        s_bar = -0.5*(sigma_s**2)
        smin = -3*sigma_s + s_bar
        smax = 3*sigma_s + s_bar
        ds = (smax - smin)/1000
        value = Z
        for y in range(0, len(n_H_mean_arr)):
            n_H_mean = n_H_mean_arr[y]
            s = np.zeros(1000)
            pdf = np.zeros(1000)
            lambda_jeans = np.zeros(1000)
            X_H2 = np.zeros(1000)
            n_H = np.zeros(1000)
            integral1 = 0
            X_H2_bar = np.zeros(1000)
            for i in range(0, 1000):
                s[i] = smin + i*ds
                n_H[i] = n_H_mean * np.exp(s[i])
                lambda_jeans[i] = make_lambda_jeans(n_H[i], s[i])
                X_H2[i] = make_X_H2(n_H[i], Z, G_o)
                pdf[i] = make_pdf(s[i], s_bar, sigma_s)
            integral1 = make_integral1()
            X_H2_bar = integral1 * X_H2
            plt.scatter(np.log10(n_H_mean), X_H2_bar[i], color=color)
    plt.xlabel('log(n_H_mean)')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H_mean) vs X_H2_bar - varying Z')
    ax.legend(  custom_lines,
                [   label + '= 0.001',
                    label + '= 0.010',
                    label + '= 1.00'  ],
                loc = 'lower right'
                    )
    plt.savefig(os.path.join('log(n_H_mean)vsX_H2_bar--Z.png'.format()))
    plt.clf()
    return s, smin, smax, sigma_s, n_H, lambda_jeans, X_H2, pdf, integral1, X_H2_bar

def varying_G_o():
    fig, ax = plt.subplots()
    mach_no = 5
    Z = 1
    G_o_arr = np.array([1., 10., 100.])
    n_H_mean_arr = np.array([1e1, 1e2, 2.5e2, 5e2, 7.5e2, 1e3, 2.5e3, 5e3, 1e4, 1e5, 1e6, 1e7])
    label = "G_o "
    color_arr = ['r','g','b']
    custom_lines = [Line2D([0], [0], color='r', lw=4),
                    Line2D([0], [0], color='g', lw=4),
                    Line2D([0], [0], color='b', lw=4)]
    for g in range(0, len(G_o_arr)):
        G_o = G_o_arr[g]
        color = str(color_arr[g])
        sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
        s_bar = -0.5*(sigma_s**2)
        smin = -3*sigma_s + s_bar
        smax = 3*sigma_s + s_bar
        ds = (smax - smin)/1000
        value = G_o
        for y in range(0, len(n_H_mean_arr)):
            n_H_mean = n_H_mean_arr[y]
            s = np.zeros(1000)
            pdf = np.zeros(1000)
            lambda_jeans = np.zeros(1000)
            X_H2 = np.zeros(1000)
            n_H = np.zeros(1000)
            integral1 = 0
            X_H2_bar = np.zeros(1000)
            for i in range(0, 1000):
                s[i] = smin + i*ds
                n_H[i] = n_H_mean * np.exp(s[i])
                lambda_jeans[i] = make_lambda_jeans(n_H[i], s[i])
                X_H2[i] = make_X_H2(n_H[i], Z, G_o)
                pdf[i] = make_pdf(s[i], s_bar, sigma_s)
            integral1 = make_integral1()
            X_H2_bar = integral1 * X_H2
            plt.scatter(np.log10(n_H_mean), X_H2_bar[i], color=color)
    plt.xlabel('log(n_H_mean)')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H_mean) vs X_H2_bar - varying G_o')
    ax.legend(  custom_lines,
                [   label + '= 1',
                    label + '= 10',
                    label + '= 100'  ],
                loc = 'lower right'
                    )
    plt.savefig(os.path.join('log(n_H_mean)vsX_H2_bar--G_o.png'.format()))
    plt.clf()
    return s, smin, smax, sigma_s, n_H, lambda_jeans, X_H2, pdf, integral1, X_H2_bar

if __name__=='__main__':
    path = 'for X_H2'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)

    # order of variables:
    # s, smin, smax, sigma_s, n_H, lambda_jeans, X_H2, pdf, integral1, X_H2_bar
    g1, g2, g3, g4, g5, g6, g7, g8, g9, g10 = varying_G_o()  #varying G_o
    m1, m2, m3, m4, m5, m6, m7, m8, m9, m10 = varying_M()  #varying mach_no
    z1, z2, z3, z4, z5, z6, z7, z8, z9, z10 = varying_Z()  #varying matallicity
