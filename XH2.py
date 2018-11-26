import numpy as np
import matplotlib.pyplot as plt
import os

def make_pdf(s):
    mach_no = 5
    n_H_mean = 100
    sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
    s_mean = -0.5*(sigma_s**2)
    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_mean)/sigma_s)**2)))
    return pdf

def make_X_H2(n_H):
    m_p = 1.672621777e-24  # g
    T_mean = 10.    #K
    K_b = 1.38064852e-16    # ergs K-1
    G = 6.67408e-8          # dyne cm^2 g^-2
    kappa = 1000 * m_p
    DC = 1.7e-11
    CC = 2.5e-17    #cm3 s-1
    rad_field_outside = 1.0     #in solar units
    Z = 1
    exp_tau = np.exp(-kappa * n_H * ((np.sqrt(K_b * T_mean / m_p)) / np.sqrt(4* np.pi * G * n_H * m_p)))
    numerator = DC * rad_field_outside * exp_tau
    denominator = CC * Z * n_H
    X_H2 = 1 / (2 + (numerator/denominator) )
    return X_H2

def plotting1(n_H, X_H2):
    path = 'for X_H2'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    #plotting X_H2 vs log(n_H)
    plt.plot(np.log10(n_H), X_H2)
    plt.xlabel('log(n_H)')
    plt.ylabel('X_H2')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs X_H2')
    plt.savefig('log(n_H)vsX_H2.png'.format())
    plt.clf()

def plotting2(n_H_bar, X_H2_bar):
    #plotting log(n_H_bar) vs X_H2_bar
    plt.plot(np.log10(n_H_bar), X_H2_bar)
    plt.xlabel('log(n_H_bar)')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H_bar) vs X_H2_bar')
    plt.savefig(os.path.join('log(n_H_bar)vsX_H2_bar.png'.format()))
    plt.clf()

if __name__=='__main__':
    exp_s = np.zeros(1000)
    pdf = np.zeros(1000)
    X_H2 = np.zeros(1000)
    n_H2 = np.zeros(1000)
    n_H_bar = np.zeros(1000)
    tot_n_H_bar = np.zeros(1000)
    tot_n_H2_bar = np.zeros(1000)
    X_H2_bar = np.zeros(1000)
    integral1 = 0
    smin = -4
    smax = 5
    s = 0
    ds = (smax - smin)/1000

    for i in range(1, 1000):
        n_H = (np.logspace(-4, 5, 1000) * i) # [H] cm^-3
        X_H2 = make_X_H2(n_H)
        n_H2 = X_H2 * n_H
        s = smin + i*ds
        pdf = make_pdf(s)
        n_H_bar += n_H * np.exp(-s)
        tot_n_H_bar += np.exp(s) * pdf * ds
        X_H2_bar += np.exp(s) * pdf * ds * X_H2
    plotting1(n_H, X_H2)
    plotting2(n_H_bar, X_H2_bar)
