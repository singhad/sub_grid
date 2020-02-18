import timing
import numpy as np
import matplotlib.pyplot as plt
import os

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

def calc_n_LW(n_H, G_o, lambda_jeans, m_p, Z):
    kappa = 1000 * m_p * Z
    rad_field_outside = G_o #in solar units
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    n_LW = rad_field_outside * exp_tau
    return n_LW

def calc_n_LW_ss(n_H, n_H2, G_o, lambda_jeans, m_p, Z):
    kappa = 1000 * m_p * Z
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
    n_LW_ss = np.zeros(1000)
    S_H2_ss = np.zeros(1000)
    N_H2_ss = np.zeros(1000)
    X_H2_ss = np.zeros(1000)
    n_H2 = np.zeros(1000)
    n_H2_ss = np.zeros(1000)
    ctr = 15
    i = 0
    n_LW = calc_n_LW(n_H, G_o, lambda_jeans, m_p, Z)
    X_H2 = calc_X_H2(n_H, Z, n_LW)
    n_H2 = n_H * X_H2
    n_H2_ss = n_H2
    while i<ctr:
        n_LW_ss, S_H2_ss, N_H2_ss = calc_n_LW_ss(n_H, n_H2, G_o, lambda_jeans, m_p, Z)
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

def calc_n_CO(n_H, X_CO, Z):
    abundance_Ctot = 1e-4 # n_C/n_H as defined by nucleosynthesis
    return n_H * abundance_Ctot * Z * X_CO # CO/cc

def all_the_stuff(M, n_H_mean, Z, G_o, count):
    T = 10 #K
    m_p = 1.672621777e-24   # g
    T_mean = 10.            #K
    K_b = 1.38064852e-16    # ergs K-1
    G = 6.67408e-8          # dyne cm^2 g^-2
    c_s = np.sqrt(K_b * T / m_p)
    s = np.zeros((10, 1000))
    pdf = np.zeros((10, 1000))
    lambda_jeans = np.zeros((10, 1000))
    n_H = np.zeros((10, 1000))
    X_CO = np.zeros((10, 1000))
    n_CO = np.zeros((10, 1000))
    X_H2 = np.zeros((10, 1000))
    n_LW = np.zeros((10, 1000))
    n_LW_ss = np.zeros((10, 1000))
    S_H2_ss = np.zeros((10, 1000))
    N_H2_ss = np.zeros((10, 1000))
    X_H2_ss = np.zeros((10, 1000))
    n_H2 = np.zeros((10, 1000))
    n_H2_ss = np.zeros((10, 1000))
    integral1 = 0
    sigma_s = np.sqrt(np.log(1 + ((0.3 * M)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -7*sigma_s + s_bar
    smax = 7*sigma_s + s_bar
    ds = (smax - smin)/1000
    for i in range(1000):
        s[count][i] = smin + i*ds
        pdf[count][i] = make_pdf(s[count][i], s_bar, sigma_s)
        n_H[count][i] = n_H_mean * np.exp(s[count][i])
        lambda_jeans[count][i] = calc_lambda_jeans(n_H[count][i], c_s, G, K_b, m_p)
        n_LW[count][i], n_LW_ss[count][i], S_H2_ss[count][i], N_H2_ss[count][i], X_H2_ss[count][i], n_H2_ss[count][i], n_H2[count][i] = self_shielding_iterations(n_H[count][i], G_o, lambda_jeans[count][i], Z, m_p)
        X_CO[count][i] = calc_X_CO(n_H[count][i], n_H2[count][i], n_LW[count][i])
        n_CO[count][i] = calc_n_CO(n_H[count][i], X_CO[count][i], Z)
    plotting(Z, M, count, n_H[count][:], pdf[count][:], lambda_jeans[count][:], X_H2_ss[count][:], X_CO[count][:], n_CO[count][:], S_H2_ss[count][:], N_H2_ss[count][:])
    #return s, pdf, n_H, lambda_jeans, n_LW, n_LW_ss, S_H2_ss, N_H2_ss, X_H2_ss, n_H2_ss, X_CO, n_CO

def plotting(Z, M, count, n_H, pdf, lambda_jeans, X_H2, X_CO, n_CO, S_H2, N_H2):
    cle = str(count)
    mcu = str(M)
    dc = str(Z)
    plt.plot(np.log10(n_H), pdf)
    plt.xlabel('log(n_H)')
    plt.ylabel('pdf')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs pdf: M=' + mcu + '; Z=' + dc)
    plt.savefig(os.path.join('log(n_H)vspdf_' + cle +'.png'.format()))
    plt.clf()

    plt.plot(np.log10(n_H), np.log10(pdf))
    plt.xlabel('log(n_H)')
    plt.ylabel('log(pdf)')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs log(pdf): M=' + mcu + '; Z=' + dc)
    plt.savefig(os.path.join('log(n_H)vslog(pdf)_' + cle +'.png'.format()))
    plt.clf()

    plt.plot(np.log10(lambda_jeans), np.log10(n_H))
    plt.xlabel('log(lambda_jeans)')
    plt.ylabel('log(n_H)')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(lambda_jeans) vs log(n_H): M=' + mcu + '; Z=' + dc)
    plt.savefig(os.path.join('log(lambda_jeans)vslog(n_H)_' + cle +'.png'.format()))
    plt.clf()

    plt.plot(np.log10(n_H), X_H2)
    plt.xlabel('log(n_H)')
    plt.ylabel('X_H2')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs X_H2: M=' + mcu + '; Z=' + dc)
    plt.savefig('log(n_H)vsX_H2_' + cle +'.png'.format())
    plt.clf()

    plt.plot(np.log10(n_CO), X_CO)
    plt.xlabel('log(n_CO)')
    plt.ylabel('X_CO')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_CO) vs X_CO: M=' + mcu + '; Z=' + dc)
    plt.savefig('log(n_CO)vsX_CO_' + cle +'.png'.format())
    plt.clf()

    plt.plot(np.log10(n_H), X_CO)
    plt.xlabel('log(n_H)')
    plt.ylabel('X_CO')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs X_CO: M=' + mcu + '; Z=' + dc)
    plt.savefig('log(n_H)vsX_CO_' + cle +'.png'.format())
    plt.clf()

    plt.plot(np.log10(N_H2), S_H2)
    plt.xlabel('log(N_H2)')
    plt.ylabel('S_H2')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(N_H2) vs S_H2: M=' + mcu + '; Z=' + dc)
    plt.savefig('log(N_H2)vsS_H2_' + cle +'.png'.format())
    plt.clf()

if __name__=='__main__':
    path = 'for X_H2_self_shielding'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    Z_arr = np.logspace(-4, -1, 10)/0.02
    mach_no_arr = np.logspace(-2, 3, 10)
    n_H_mean_arr = np.logspace(-5, 4, 10)
    for m in range(0, len(mach_no_arr)):
        mach_no = mach_no_arr[m]
        n_H_mean = n_H_mean_arr[m]
        Z = Z_arr[m]
        G_o = 1
        #s, pdf, n_H, lambda_jeans, n_LW, n_LW_ss, S_H2_ss, N_H2_ss, X_H2_ss, n_H2_ss, X_CO, n_CO =
        all_the_stuff(mach_no, n_H_mean, Z, G_o, m)
