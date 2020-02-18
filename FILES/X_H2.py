import numpy as np
import matplotlib.pyplot as plt
import os

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
    smin = -4*sigma_s + s_bar
    smax = 4*sigma_s + s_bar
    ds = (smax - smin)/1000
    n_H_mean = 1e4
    Z = 1
    G_o = 1
    for i in range(0, 1000):
        s[i] = smin + i*ds
        n_H = n_H_mean * np.exp(s)
        lambda_jeans = make_lambda_jeans(n_H, s)
        X_H2 = make_X_H2(n_H, Z, G_o)
        pdf = make_pdf(s, s_bar, sigma_s)
        integral1 += np.exp(s[i]) * pdf[i] * ds   #this should be ~1
    plotting(n_H, pdf, lambda_jeans, X_H2)
    return n_H, pdf, lambda_jeans, X_H2, integral1

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

if __name__=='__main__':
    path = 'for H2_density'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    n_H, pdf, lambda_jeans, X_H2, integral1 = make_integral1()
