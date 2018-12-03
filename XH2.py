import numpy as np
import matplotlib.pyplot as plt
import os

def make_pdf(s, s_bar, sigma_s):
    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_bar)/sigma_s)**2)))
    return pdf

def make_lambda_jeans(n_H, s):
    m_p = 1.672621777e-24  # g
    T_mean = 10.    #K
    K_b = 1.38064852e-16    # ergs K-1
    G = 6.67408e-8          # dyne cm^2 g^-2
    lambda_jeans = ((np.sqrt(K_b * T_mean / m_p)) / np.sqrt(4* np.pi * G * n_H * m_p))
    return lambda_jeans

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

def plotting1(s, pdf):
    plt.plot(s, pdf)
    plt.xlabel('s')
    plt.ylabel('pdf')
    plt.grid(b=True, which='both', axis='both')
    plt.title('s vs pdf')
    plt.savefig(os.path.join('Svspdf.png'.format()))
    plt.clf()
    plt.plot(s, np.log10(pdf))
    plt.xlabel('s')
    plt.ylabel('log(pdf)')
    plt.grid(b=True, which='both', axis='both')
    plt.title('s vs log(pdf)')
    plt.savefig(os.path.join('Svslog(pdf).png'.format()))
    plt.clf()

def plotting2(s, lambda_jeans):
    plt.plot(np.log10(lambda_jeans), s)
    plt.xlabel('log(lambda_jeans)')
    plt.ylabel('s')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(lambda_jeans) vs s')
    plt.savefig(os.path.join('log(lambda_jeans)vsS.png'.format()))
    plt.clf()

def plotting3(s, X_H2):
    plt.plot(s, X_H2)
    plt.xlabel('s')
    plt.ylabel('X_H2')
    plt.grid(b=True, which='both', axis='both')
    plt.title('s vs X_H2')
    plt.savefig('SvsX_H2.png'.format())
    plt.clf()

def plotting4(s, X_H2_bar):
    plt.plot(s, X_H2_bar)
    plt.xlabel('s')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('s vs X_H2_bar')
    plt.savefig(os.path.join('SvsX_H2_bar.png'.format()))
    plt.clf()

def plotting5(s_bar, X_H2_bar):
    plt.scatter(s_bar, X_H2_bar)
    plt.xlabel('s_bar')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('s_bar vs X_H2_bar')

def plotting6(n_H_mean, X_H2_bar):
    plt.scatter(np.log10(n_H_mean), X_H2_bar)
    plt.xlabel('log(n_H_mean)')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H_mean) vs X_H2_bar')

def function1():
    s = np.zeros(1000)
    pdf = np.zeros(1000)
    mach_no = 5
    sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -3*sigma_s + s_bar
    smax = 3*sigma_s + s_bar
    ds = (smax - smin)/1000
    n_H_mean = 1e4
    for i in range(0, 1000):
        s[i] = smin + i*ds
        pdf = make_pdf(s, s_bar, sigma_s)
    plotting1(s, pdf)

def function23():
    s = np.zeros(1000)
    X_H2 = np.zeros(1000)
    n_H = np.zeros(1000)
    lambda_jeans = np.zeros(1000)
    mach_no = 5
    sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -3*sigma_s + s_bar
    smax = 3*sigma_s + s_bar
    ds = (smax - smin)/1000
    n_H_mean = 1e4
    for i in range(0, 1000):
        s[i] = smin + i*ds
        n_H = n_H_mean * np.exp(s)
        lambda_jeans = make_lambda_jeans(n_H, s)
        X_H2 = make_X_H2(n_H)
    plotting2(s, lambda_jeans)
    plotting3(s, X_H2)
    return smin, smax, s, lambda_jeans, n_H

def function4():
    s = np.zeros(1000)
    pdf = np.zeros(1000)
    lambda_jeans = np.zeros(1000)
    X_H2 = np.zeros(1000)
    n_H = np.zeros(1000)
    tot_n_H_bar = np.zeros(1000)
    X_H2_bar = np.zeros(1000)
    mach_no = 5
    sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -3*sigma_s + s_bar
    smax = 3*sigma_s + s_bar
    ds = (smax - smin)/1000
    n_H_mean = 1e4
    s[0] = smin
    for i in range(0, 1000):
        s[i] = smin + i*ds
        n_H = n_H_mean * np.exp(s)
        X_H2 = make_X_H2(n_H)
        pdf = make_pdf(s, s_bar, sigma_s)
        tot_n_H_bar += np.exp(s) * pdf * ds
        X_H2_bar += np.exp(s) * pdf * ds * X_H2
    plotting4(s, X_H2_bar)
    return smin, smax, s, n_H, tot_n_H_bar, X_H2_bar

def function5():
    s = np.zeros(1000)
    pdf = np.zeros(1000)
    X_H2 = np.zeros(1000)
    n_H = np.zeros(1000)
    tot_n_H_bar = np.zeros(1000)
    X_H2_bar = np.zeros(1000)
    mach_no_arr = np.array([1., 5., 10., 15., 20., 50., 100.])
    for m in range(0, len(mach_no_arr)):
        mach_no = mach_no_arr[m]
        sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
        s_bar = -0.5*(sigma_s**2)
        smin = -3*sigma_s + s_bar
        smax = 3*sigma_s + s_bar
        ds = (smax - smin)/1000
        n_H_mean = 1e4
        s[0] = smin
        #s_mean = [s_bar]*len(X_H2_bar)
        for i in range(0, 1000):
            s[i] = smin + i*ds
            n_H[i] = n_H_mean * np.exp(s[i])
            X_H2[i] = make_X_H2(n_H[i])
            pdf[i] = make_pdf(s[i], s_bar, sigma_s)
            tot_n_H_bar[i] = np.exp(s[i]) * pdf[i] * ds
            X_H2_bar[i] = 0.9996868549 * X_H2[i]
        plotting5(s_bar, X_H2_bar[i])
    plt.savefig(os.path.join('S_barvsX_H2_bar--5.png'.format()))
    plt.clf()
    return tot_n_H_bar, X_H2_bar

def function6():
    s = np.zeros(1000)
    pdf = np.zeros(1000)
    X_H2 = np.zeros(1000)
    n_H = np.zeros(1000)
    tot_n_H_bar = np.zeros(1000)
    X_H2_bar = np.zeros(1000)
    n_H_mean_arr = np.array([1e1, 1e2, 2.5e2, 5e2, 7.5e2, 1e3, 2.5e3, 5e3, 1e4, 1e5, 1e6, 1e7])
    for m in range(0, len(n_H_mean_arr)):
        n_H_mean = n_H_mean_arr[m]
        mach_no = 5
        sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
        s_bar = -0.5*(sigma_s**2)
        smin = -3*sigma_s + s_bar
        smax = 3*sigma_s + s_bar
        ds = (smax - smin)/1000
        #s_mean = [s_bar]*len(X_H2_bar)
        for i in range(0, 1000):
            s[i] = smin + i*ds
            n_H[i] = n_H_mean * np.exp(s[i])
            X_H2[i] = make_X_H2(n_H[i])
            pdf[i] = make_pdf(s[i], s_bar, sigma_s)
            tot_n_H_bar[i] = np.exp(s[i]) * pdf[i] * ds
            X_H2_bar[i] = 0.9996868549 * X_H2[i]
        plotting6(n_H_mean, X_H2_bar[i])
    plt.savefig(os.path.join('log(n_H_mean)vsX_H2_bar--6.png'.format()))
    plt.clf()
    return tot_n_H_bar, X_H2_bar


if __name__=='__main__':
    path = 'for X_H2'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)

    #function1() #SvsPDF
    #smin, smax, s, lambda_jeans, n_H = function23()  #SvsX_H2 & #log(lambda_jeans)vsS
    smin, smax, s, n_H, tot_n_H_bar, X_H2_bar = function4()  #SvsX_H2_bar
    #tot_n_H_bar, X_H2_bar = function5()  #S_barvsX_H2_bar : varying mach_no
    #tot_n_H_bar, X_H2_bar = function6()  #n_H_meanvsX_H2_bar : varying n_H_mean
