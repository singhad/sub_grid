import numpy as np
import matplotlib.pyplot as plt
import os

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

def calc_n_LW_ss(n_H, n_H2, G_o, lambda_jeans):
    m_p = 1.672621777e-24   # g
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

def CO_density():
    n_LW_1 = np.zeros(1000)
    n_LW_2 = np.zeros(1000)
    n_LW_3 = np.zeros(1000)
    n_LW_4 = np.zeros(1000)
    n_LW_5 = np.zeros(1000)
    n_LW_6 = np.zeros(1000)
    n_LW_7 = np.zeros(1000)
    n_LW_8 = np.zeros(1000)
    n_LW_9 = np.zeros(1000)
    n_LW_10 = np.zeros(1000)

    X_H2_1 = np.zeros(1000)
    X_H2_2 = np.zeros(1000)
    X_H2_3 = np.zeros(1000)
    X_H2_4 = np.zeros(1000)
    X_H2_5 = np.zeros(1000)
    X_H2_6 = np.zeros(1000)
    X_H2_7 = np.zeros(1000)
    X_H2_8 = np.zeros(1000)
    X_H2_9 = np.zeros(1000)
    X_H2_10 = np.zeros(1000)

    n_H2_1 = np.zeros(1000)
    n_H2_2 = np.zeros(1000)
    n_H2_3 = np.zeros(1000)
    n_H2_4 = np.zeros(1000)
    n_H2_5 = np.zeros(1000)
    n_H2_6 = np.zeros(1000)
    n_H2_7 = np.zeros(1000)
    n_H2_8 = np.zeros(1000)
    n_H2_9 = np.zeros(1000)
    n_H2_10 = np.zeros(1000)

    S_H2_1 = np.zeros(1000)
    S_H2_2 = np.zeros(1000)
    S_H2_3 = np.zeros(1000)
    S_H2_4 = np.zeros(1000)
    S_H2_5 = np.zeros(1000)
    S_H2_6 = np.zeros(1000)
    S_H2_7 = np.zeros(1000)
    S_H2_8 = np.zeros(1000)
    S_H2_9 = np.zeros(1000)
    S_H2_10 = np.zeros(1000)

    N_H2_1 = np.zeros(1000)
    N_H2_2 = np.zeros(1000)
    N_H2_3 = np.zeros(1000)
    N_H2_4 = np.zeros(1000)
    N_H2_5 = np.zeros(1000)
    N_H2_6 = np.zeros(1000)
    N_H2_7 = np.zeros(1000)
    N_H2_8 = np.zeros(1000)
    N_H2_9 = np.zeros(1000)
    N_H2_10 = np.zeros(1000)

    s = np.zeros(1000)
    pdf = np.zeros(1000)
    lambda_jeans = np.zeros(1000)
    X_H2 = np.zeros(1000)
    n_H = np.zeros(1000)
    n_LW = np.zeros(1000)
    n_LW_ss = np.zeros(1000)
    S_H2 = np.zeros(1000)
    n_H2 = np.zeros(1000)
    X_CO = np.zeros(1000)
    n_CO = np.zeros(1000)
    N_H2 = np.zeros(1000)
    mach_no = 5
    sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -4*sigma_s + s_bar
    smax = 4*sigma_s + s_bar
    ds = (smax - smin)/1000
    n_H_mean = 1e4
    Z = 1
    G_o = 1
    X_H2_a = np.zeros(1000)
    n_H2_a = np.zeros(1000)

    for i in range(0, 1000):
        s[i] = smin + i*ds

    pdf = make_pdf(s, s_bar, sigma_s)
    n_H = n_H_mean * np.exp(s)
    lambda_jeans = calc_lambda_jeans(n_H)

    n_LW = calc_n_LW(n_H, G_o, lambda_jeans)
    X_H2_a = calc_X_H2(n_H, Z, n_LW)
    n_H2_a = n_H * X_H2_a
    a_1 = X_H2_a
    a_2 = n_H2_a

    n_LW_1, S_H2_1, N_H2_1 = calc_n_LW_ss(n_H, n_H2_a, G_o, lambda_jeans)
    X_H2_1 = calc_X_H2(n_H, Z, n_LW_1)
    n_H2_1 = n_H * X_H2_1

    n_LW_2, S_H2_2, N_H2_2 = calc_n_LW_ss(n_H, n_H2_1, G_o, lambda_jeans)
    X_H2_2 = calc_X_H2(n_H, Z, n_LW_2)
    n_H2_2 = n_H * X_H2_2

    n_LW_3, S_H2_3, N_H2_3 = calc_n_LW_ss(n_H, n_H2_2, G_o, lambda_jeans)
    X_H2_3 = calc_X_H2(n_H, Z, n_LW_3)
    n_H2_3 = n_H * X_H2_3

    n_LW_4, S_H2_4, N_H2_4 = calc_n_LW_ss(n_H, n_H2_3, G_o, lambda_jeans)
    X_H2_4 = calc_X_H2(n_H, Z, n_LW_4)
    n_H2_4 = n_H * X_H2_4

    n_LW_5, S_H2_5, N_H2_5 = calc_n_LW_ss(n_H, n_H2_4, G_o, lambda_jeans)
    X_H2_5 = calc_X_H2(n_H, Z, n_LW_5)
    n_H2_5 = n_H * X_H2_5

    n_LW_6, S_H2_6, N_H2_6 = calc_n_LW_ss(n_H, n_H2_5, G_o, lambda_jeans)
    X_H2_6 = calc_X_H2(n_H, Z, n_LW_6)
    n_H2_6 = n_H * X_H2_6

    n_LW_7, S_H2_7, N_H2_7 = calc_n_LW_ss(n_H, n_H2_6, G_o, lambda_jeans)
    X_H2_7 = calc_X_H2(n_H, Z, n_LW_7)
    n_H2_7 = n_H * X_H2_7

    n_LW_8, S_H2_8, N_H2_8 = calc_n_LW_ss(n_H, n_H2_7, G_o, lambda_jeans)
    X_H2_8 = calc_X_H2(n_H, Z, n_LW_8)
    n_H2_8 = n_H * X_H2_8

    n_LW_9, S_H2_9, N_H2_9 = calc_n_LW_ss(n_H, n_H2_8, G_o, lambda_jeans)
    X_H2_9 = calc_X_H2(n_H, Z, n_LW_9)
    n_H2_9 = n_H * X_H2_9

    n_LW_10, S_H2_10, N_H2_10 = calc_n_LW_ss(n_H, n_H2_9, G_o, lambda_jeans)
    X_H2_10 = calc_X_H2(n_H, Z, n_LW_10)
    n_H2_10 = n_H * X_H2_10

    n_LW_ss, S_H2, N_H2 = calc_n_LW_ss(n_H, n_H2_10, G_o, lambda_jeans)
    X_H2 = calc_X_H2(n_H, Z, n_LW_ss)
    n_H2 = n_H * X_H2

    X_CO = calc_X_CO(n_H, n_H2_a, n_LW)
    n_CO = calc_n_CO(n_H, X_CO)

    plotting(n_H, pdf, lambda_jeans, X_H2, X_CO, n_CO, S_H2, N_H2)
    return n_H, n_H2, X_H2, n_LW, n_LW_ss, S_H2, N_H2, n_LW_1, n_LW_2, n_LW_3, n_LW_4, n_LW_5, n_LW_6, n_LW_7, n_LW_8, n_LW_9, n_LW_10, X_H2_1, X_H2_2, X_H2_3, X_H2_4, X_H2_5, X_H2_6, X_H2_7, X_H2_8, X_H2_9, X_H2_10, n_H2_1, n_H2_2, n_H2_3, n_H2_4, n_H2_5, n_H2_6, n_H2_7, n_H2_8, n_H2_9, n_H2_10, X_H2_a, n_H2_a

def plotting(n_H, pdf, lambda_jeans, X_H2, X_CO, n_CO, S_H2, N_H2):
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

    plt.plot(np.log10(n_CO), X_CO)
    plt.xlabel('log(n_CO)')
    plt.ylabel('X_CO')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_CO) vs X_CO')
    plt.savefig('log(n_CO)vsX_CO.png'.format())
    plt.clf()

    plt.plot(np.log10(n_H), X_CO)
    plt.xlabel('log(n_H)')
    plt.ylabel('X_CO')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs X_CO')
    plt.savefig('log(n_H)vsX_CO.png'.format())
    plt.clf()

    plt.plot(np.log10(N_H2), S_H2)
    plt.xlabel('log(N_H2)')
    plt.ylabel('S_H2')
    plt.grid(b=True, which='both', axis='both')
    plt.title("log(N_H2) vs S_H2")
    plt.savefig('log(N_H2)vsS_H2.png'.format())
    plt.clf()

if __name__=='__main__':
    path = 'for X_H2_self_shielding'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    n_H, n_H2, X_H2, n_LW, n_LW_ss, S_H2, N_H2, n_LW_1, n_LW_2, n_LW_3, n_LW_4, n_LW_5, n_LW_6, n_LW_7, n_LW_8, n_LW_9, n_LW_10, X_H2_1, X_H2_2, X_H2_3, X_H2_4, X_H2_5, X_H2_6, X_H2_7, X_H2_8, X_H2_9, X_H2_10, n_H2_1, n_H2_2, n_H2_3, n_H2_4, n_H2_5, n_H2_6, n_H2_7, n_H2_8, n_H2_9, n_H2_10, X_H2_a, n_H2_a = CO_density()
