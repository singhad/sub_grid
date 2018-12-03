'''
    for a in range(0, len(n_H_mean_arr)):
        n_H_mean = n_H_mean_arr[a]
        for i in range(0, 1000):
            s[i] = smin + i*ds
            n_H[i] = n_H_mean * np.exp(s[i])
            X_H2[i] = make_X_H2(n_H[i])
            pdf[i] = make_pdf(s[i], s_bar, sigma_s)
            tot_n_H_bar[i] = np.exp(s[i]) * pdf[i] * ds
            X_H2_bar[i] = np.exp(s[i]) * pdf[i] * ds * X_H2[i]
        #X_H2_bar_arr[a] = X_H2_bar
    plotting3(s, pdf)

'''

"""
# ---------------
# Code to find integration and also plot multiple plots
# ---------------

import numpy as np
import matplotlib.pyplot as plt
import os

def make_pdf(s, sigma_s, s_mean):
    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_mean)/sigma_s)**2)))
    return pdf

def make_X_H2(n_H, Z):
    m_p = 1.672621777e-24  # g
    T_mean = 10.    #K
    K_b = 1.38064852e-16    # ergs K-1
    G = 6.67408e-8          # dyne cm^2 g^-2
    kappa = 1000 * m_p
    DC = 1.7e-11
    CC = 2.5e-17    #cm3 s-1
    rad_field_outside = 1.0     #in solar units
    exp_tau = np.exp(-kappa * n_H * ((np.sqrt(K_b * T_mean / m_p)) / np.sqrt(4* np.pi * G * n_H * m_p)))
    numerator = DC * rad_field_outside * exp_tau
    denominator = CC * Z * n_H
    X_H2 = 1 / (2 + (numerator/denominator) )
    return X_H2

def plot_X_H2_bar_vs_n_H_bar(n_H_bar, X_H2_bar, label, value):
    plt.plot(np.log10(n_H_bar), X_H2_bar, label =(label + '= ' + str(value)))
    plt.xlabel('log(n_H_bar)')
    plt.ylabel('X_H2_bar')
    plt.legend()
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H_bar) vs X_H2_bar - varying ' + label)

if __name__=='__main__':
    path = 'multiple_plots'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    mach_no_arr = np.array([1., 5., 10., 15., 20., 50., 100.])
    Z_arr = np.array([0.001, 0.01, 0.1, 1., 10., 100.])
    n_H_mean_arr = np.array([1e-1, 1, 1e1, 1e2, 1e3, 1e4])
    pdf = np.zeros(1000)
    X_H2 = np.zeros(1000)
    n_H2 = np.zeros(1000)
    n_H_bar = np.zeros(1000)
    tot_n_H_bar = np.zeros(1000)
    X_H2_bar = np.zeros(1000)
    smin = -4
    smax = 5
    ds = (smax - smin)/1000

    choice = [1, 2, 3]
    for ch in choice:
        value = 0
        label = ""
        if ch == 1:
            mach_no = 5
            sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
            s_mean = -0.5*(sigma_s**2)
            n_H_mean = 100
            label = "Metallicity "
            for a in range(0, len(Z_arr)):
                Z = Z_arr[a]
                value = Z
                tot_n_H_bar = np.zeros(1000)
                X_H2 = np.zeros(1000)
                for i in range(1, 1000):
                    n_H = (np.logspace(-4, 5, 1000) * i)
                    X_H2 = make_X_H2(n_H, Z)
                    n_H2 = X_H2 * n_H
                    s = smin + i*ds
                    n_H_bar = n_H * np.exp(-s)
                    pdf = make_pdf(s, sigma_s, s_mean)
                    tot_n_H_bar += np.exp(s) * pdf * ds
                    X_H2_bar = tot_n_H_bar * X_H2
                plot_X_H2_bar_vs_n_H_bar(n_H_bar, X_H2_bar, label, value)
            plt.savefig('log(n_H_bar)vsX_H2_bar__Z.png'.format())
            plt.clf()
        elif ch == 2:
            Z = 1
            n_H_mean = 100
            label = "Mach no "
            for b in range(0, len(mach_no_arr)):
                mach_no = mach_no_arr[b]
                sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
                s_mean = -0.5*(sigma_s**2)
                value = mach_no
                tot_n_H_bar = np.zeros(1000)
                X_H2 = np.zeros(1000)
                for i in range(1, 1000):
                    n_H = (np.logspace(-4, 5, 1000) * i)
                    X_H2 = make_X_H2(n_H, Z)
                    n_H2 = X_H2 * n_H
                    s = smin + i*ds
                    n_H_bar = n_H * np.exp(-s)
                    pdf = make_pdf(s, sigma_s, s_mean)
                    tot_n_H_bar += np.exp(s) * pdf * ds
                    X_H2_bar = tot_n_H_bar * X_H2
                plot_X_H2_bar_vs_n_H_bar(n_H_bar, X_H2_bar, label, value)
            plt.savefig('log(n_H_bar)vsX_H2_bar__M.png'.format())
            plt.clf()
        elif ch == 3:
            mach_no = 5
            sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
            s_mean = -0.5*(sigma_s**2)
            Z = 1
            label = "n_H_mean "
            for c in range(0, len(n_H_mean_arr)):
                n_H_mean = n_H_mean_arr[c]
                value = n_H_mean
                tot_n_H_bar = np.zeros(1000)
                X_H2 = np.zeros(1000)
                for i in range(1, 1000):
                    n_H = (np.logspace(-4, 5, 1000) * i)
                    X_H2 = make_X_H2(n_H, Z)
                    n_H2 = X_H2 * n_H
                    s = smin + i*ds
                    n_H_bar = n_H * np.exp(-s)
                    pdf = make_pdf(s, sigma_s, s_mean)
                    tot_n_H_bar += np.exp(s) * pdf * ds
                    X_H2_bar = tot_n_H_bar * X_H2
                plot_X_H2_bar_vs_n_H_bar(n_H_bar, X_H2_bar, label, value)
            plt.savefig('log(n_H_bar)vsX_H2_bar__n_H_mean.png'.format())
            plt.clf()

"""


# ---------------
# Working code for integration to find X_H2_bar
# ---------------

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
        n_H_bar = n_H * np.exp(-s)
        tot_n_H_bar += np.exp(s) * pdf * ds
        X_H2_bar += np.exp(s) * pdf * ds * X_H2
    plotting1(n_H, X_H2)
    plotting2(n_H_bar, X_H2_bar)






#*******************************************************************************

'''
# ---------------
# Code for the plotting funcion that stores multiple images
# ---------------
def plotting(x, pdf, lambda_jeans_cm, n_H, X_H2, X_H2_bar):
    path = 'with_integration'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    #plotting log(n) vs log(PDF)
    plt.plot(x, np.log(pdf))
    plt.xlabel('log(n) [H/cc]')
    plt.ylabel('log(PDF)')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n) vs log(PDF)')
    plt.savefig(os.path.join('XvslogPDF','image_%04d.png'.format()))
    #plt.savefig('XvslogPDFimage.png'.format())
    plt.clf()

    #plotting log(n) vs PDF
    plt.plot(x, pdf)
    plt.xlabel('log(n) [H/cc]')
    plt.ylabel('PDF')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n) vs PDF')
    plt.savefig(os.path.join('XvsPDF','image_%04d.png'.format()))
    #plt.savefig('XvsPDF.png'.format())
    plt.clf()

    #plotting log of Jeans length in cm log(lambda_jeans_cm) vs log(n)
    plt.plot(np.log(lambda_jeans_cm), x)
    plt.ylabel('log(n) [H/cc]')
    plt.xlabel('log(lambda_Jeans) [cm]')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(lambda_Jeans) vs log(n)')
    plt.savefig(os.path.join('log(lambda_Jeans)vsX','image_%04d.png'.format()))
    #plt.savefig('log(lambda_Jeans)vsX.png'.format())
    plt.clf()

    #plotting X_H2 vs log(n_H)
    plt.plot(x, X_H2)
    plt.xlabel('log(n_H) [H/cc]')
    plt.ylabel('X_H2')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs X_H2')
    plt.savefig(os.path.join('log(n_H)vsX_H2','image_%04d.png'.format()))
    #plt.savefig('log(n_H)vsX_H2.png'.format())
    plt.clf()

    #plotting n_H vs X_H2_bar
    plt.plot(n_H, X_H2_bar)
    plt.xlabel('n_H [H/cc]')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('n_H vs X_H2_bar')
    plt.savefig(os.path.join('n_HvsX_H2_bar','image_%04d.png'.format()))
    #plt.savefig('n_HvsX_H2_bar.png'.format())
    plt.clf()
'''


'''
# ---------------
# Code for just the for loop - with integration
# ---------------

   for n_H_mean_ctr in range(len(n_H_mean_arr)-1):
        n_H_mean = n_H_mean_arr[n_H_mean_ctr] # [H] cm^-3
        for i in range(1, n_H_range):
            n_H = (np.logspace(-4, 5, 100) * i) # [H] cm^-3
            x = np.log(n_H/n_H_mean)
            pdf = make_pdf(x, x_mean, s)
            lambda_jeans_cm = calc_lambda_jeans(T_mean, n_H)
            tau = calc_optical_depth(n_H, lambda_jeans_cm)
            n_LW = calc_num_LW(tau)
            X_H2 = calc_X_H2(n_H, n_LW, Z)
            n_H2 = calc_num_H2(n_H, X_H2)
            n_H2_bar = n_H * pdf * X_H2
            for s, item2 in enumerate(n_H2_bar):
                for k, item in enumerate(n_H2_bar):
                    tot_n_H2_bar[k] = quad(calc_n_H2_bar, np.min(n_H), np.max(n_H), args=(n_H2_bar[k], ))[0]
            X_H2_bar =  n_H_mean/n_H2_bar
            tot_X_H2_bar = calc_X_H2_bar(n_H_mean, tot_n_H2_bar)'''
'''
# ---------------
# Code attempting to create multi-dimensional array for plotting,
# also containing the plotting function for plotting all different values
# ---------------
def plotting(x, pdf, lambda_jeans_cm, n_H, X_H2, X_H2_bar):
    #plotting log(n) vs log(PDF)
    plt.plot(x, np.log(pdf), lw=1)
    plt.xlabel('log(n) [H/cc]')
    plt.ylabel('log(PDF)')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n) vs log(PDF)')
    plt.savefig('XvslogPDFimage_{k:06d}.png'.format(k=k))
    plt.clf()
    #plotting log(n) vs PDF
    plt.plot(x, pdf, lw=1)
    plt.xlabel('log(n) [H/cc]')
    plt.ylabel('PDF')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n) vs PDF')
    plt.savefig('XvsPDF_{k:06d}.png'.format(k=k))
    plt.clf()
    #plotting log of Jeans length in cm log(lambda_jeans_cm) vs log(n)
    plt.plot(np.log(lambda_jeans_cm), x, lw=1)
    plt.ylabel('log(n) [H/cc]')
    plt.xlabel('log(lambda_Jeans) [cm]')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(lambda_Jeans) vs log(n)')
    plt.savefig('log(lambda_Jeans)vsX_{k:06d}.png'.format(k=k))
    plt.clf()
    #plotting X_H2 vs log(n_H)
    plt.plot(x, X_H2, lw=1)
    plt.xlabel('log(n_H) [H/cc]')
    plt.ylabel('X_H2')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs X_H2')
    plt.savefig('log(n_H)vsX_H2_{k:06d}.png'.format(k=k))
    plt.clf()
    #plotting n_H vs X_H2_bar
    plt.plot(n_H, X_H2_bar, lw=1)
    plt.xlabel('n_H [H/cc]')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('n_H vs X_H2_bar')
    plt.savefig('n_HvsX_H2_bar_{k:06d}.png'.format(k=k))
    plt.clf()

if __name__=='__main__':
    fig, ax = plt.subplots()
    m_p = 1.672621777e-24  # g
    T_mean = 10.    #K
    mach_no_arr = np.array([1.1, 5., 10., 15.])
    Z_arr = np.array([0.0011, 0.011, 0.11, 1.1])
    n_H_mean_arr = np.array([1e-3, 1e-2, 1e-1, 1e1, 1e2, 1e3])
    n_H_range = 100
    x_mean = 1
    n_H = np.zeros([4, 4, 6, 100])
    pdf = np.zeros([4, 4, 6, 100])
    x = np.zeros([4, 4, 6, 100])
    lambda_jeans_cm = np.zeros([4, 4, 6, 100])    # cm
    tau = np.zeros([4, 4, 6, 100])     #optical depth
    n_LW = np.zeros([4, 4, 6, 100])     #number of Lyman-Werner photons
    n_H2 = np.zeros([4, 4, 6, 100])
    X_H2 = np.zeros([4, 4, 6, 100])
    n_H2_bar = np.zeros([4, 4, 6, 100])
    X_H2_bar = np.zeros([4, 4, 6, 100])

    n_H_mean = 100 # [H] cm^-3
    for i in range(len(Z_arr)):
        Z = Z_arr[i]
        for j in range(len(mach_no_arr)):
            mach_no = mach_no_arr[j]
            s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2))) #s in pdf
            for k in range(len(n_H_mean_arr)):
                n_H_mean = n_H_mean_arr[k] # [H] cm^-3
                for l in range(1, n_H_range):
                    n_H[i][j][k] = (np.logspace(-3, 3, 100) * l) # [H] cm^-3
                    x[i][j][k][l] = np.log(n_H[i][j][k][l]/n_H_mean)
                    pdf[i][j][k][l] = make_pdf(x[i][j][k][l], x_mean, s)
                    lambda_jeans_cm[i][j][k][l] = calc_lambda_jeans(T_mean, n_H[i][j][k][l])
                    tau[i][j][k][l] = calc_optical_depth(n_H[i][j][k][l], lambda_jeans_cm[i][j][k][l])
                    n_LW[i][j][k][l] = calc_num_LW(tau[i][j][k][l])
                    X_H2[i][j][k][l] = calc_X_H2(n_H[i][j][k][l], n_LW[i][j][k][l], Z)
                    n_H2[i][j][k][l] = calc_num_H2(n_H[i][j][k][l], X_H2[i][j][k][l])
                    n_H2_bar[i][j][k][l] = n_H[i][j][k][l] * pdf[i][j][k][l] * X_H2[i][j][k][l]
                    X_H2_bar[i][j][k][l] =  n_H_mean/n_H2_bar[i][j][k][l]

'''



'''
# ---------------
# Code using classes
# ---------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad


class H2_density:

    m_p = 1.672621777e-24  # g
    T_mean = 10.    #K
    mach_no_arr = np.array([1., 5., 10., 15.])
    Z_arr = np.array([0.001, 0.01, 0.1, 1])
    n_H_mean_arr = np.array([1e-3, 1e-2, 1e-1, 1e1, 1e2, 1e3])
    n_H_range = 100
    x_mean = 1
    pdf = np.zeros(100)
    x = np.zeros(100)
    lambda_jeans_cm = np.zeros(100)    # cm
    tau = np.zeros(100)     #optical depth
    n_LW = np.zeros(100)    #number of Lyman-Werner photons
    n_H2 = np.zeros(100)
    X_H2 = np.zeros(100)
    n_H_bar = np.zeros(100)
    n_H2_bar = np.zeros(100)
    n_H2_bar = np.zeros(100)
    X_H2_bar = np.zeros(100)

    self.Z_arr = Z_arr
    #function to generate the PDF
    def make_pdf(x, x_mean, s):
        pdf = (1/np.sqrt(2*np.pi*(s**2))) * (np.exp(-0.5*(((x - x_mean)/s)**2)))
        return pdf

    #function to calculate Jeans length
    def calc_lambda_jeans(T_mean, n_H):
        K_b = 1.38064852e-16    # ergs K-1
        m_p = 1.672621777e-24   # g
        G = 6.67408e-8          # dyne cm^2 g^-2
        c_s = np.sqrt(K_b * T_mean / m_p)       # cm s^-1
        lambda_jeans_cm = c_s / np.sqrt(4* np.pi * G * n_H * m_p)        # cm
        return lambda_jeans_cm

    #function to return optical depth
    def calc_optical_depth(n_H, lambda_jeans_cm):
        kappa_dust = 1000       # cm^2 g^-1
        m_p = 1.672621777e-24   # g
        kappa = kappa_dust * m_p
        tau = kappa * n_H * lambda_jeans_cm
        return tau

    def calc_num_LW(tau):
        rad_field_outside = 1.0     #in solar units
        num_LW = rad_field_outside * np.exp(-tau)
        return num_LW

    def calc_X_H2(n_H, n_LW, Z):
        destruction_coeff = 1.7e-11
        creation_coeff = 2.5e-17    #cm3 s-1
        factor = (destruction_coeff * n_LW) / (creation_coeff * Z * n_H)
        X_H2 = 1 / (2 + factor)
        return X_H2

    def calc_num_H2(n_H, X_H2):
        n_H2 = n_H * X_H2   #H2/cc
        return n_H2

    def calc_n_H_bar():
        n_H_bar = n_H * pdf
        return n_H_bar

    def calc_n_H2_bar():
        n_H2_bar = X_H2 * n_H * pdf
        return n_H2_bar

    def calc_X_H2_bar():
        X_H2_bar = n_H_bar/n_H2_bar
        return X_H2_bar

    def plotting(self):
        #variable_declaration()
        fig, ax = plt.subplots()
        for z_ctr in range(len(self.Z_arr)-1):
            Z = Z_arr[z_ctr]
            for mach_ctr in range(len(mach_no_arr)-1):
                mach_no = mach_no_arr[mach_ctr]
                s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2))) #s in pdf
                for n_H_mean_ctr in range(len(n_H_mean_arr)-1):
                    n_H_mean = n_H_mean_arr[n_H_mean_ctr] # [H] cm^-3
                    for i in range(1, n_H_range):
                        n_H = (np.logspace(-4, 5, 100) * i) # [H] cm^-3
                        x = np.log(n_H/n_H_mean)
                        pdf = make_pdf(x, x_mean, s)
                        lambda_jeans_cm = calc_lambda_jeans(T_mean, n_H)
                        tau = calc_optical_depth(n_H, lambda_jeans_cm)
                        n_LW = calc_num_LW(tau)
                        X_H2 = calc_X_H2(n_H, n_LW, Z)
                        n_H2 = calc_num_H2(n_H, X_H2)
                        #n_H_bar = quad(calc_n_H_bar, 0, 100, args=pdf, )
                        #n_H2_bar = quad(calc_n_H2_bar, 0, 100, args=(pdf, X_H2, ))
                        #X_H2_bar = calc_X_H2_bar(n_H_bar, n_H2_bar)
                    #plotting log(n) vs log(PDF)
                    plt.plot(x, np.log(pdf), lw=1, color='b')
                    plt.xlabel('log(n) [H/cc]')
                    plt.ylabel('log(PDF)')
                    plt.grid(b=True, which='both', axis='both')
                    #plt.legend()
                    plt.title('log(n) vs log(PDF)')
                    plt.savefig('XvslogPDF_fromcode.png'.format(i=i))
                    plt.clf()
                    #plotting log(n) vs PDF
                    plt.plot(x, pdf, lw=1, color='b')
                    plt.xlabel('log(n) [H/cc]')
                    plt.ylabel('PDF')
                    plt.grid(b=True, which='both', axis='both')
                    #plt.legend()
                    plt.title('log(n) vs PDF')
                    plt.savefig('XvsPDF_fromcode.png'.format(i=i))
                    plt.clf()
                    #plotting log of Jeans length in cm log(lambda_jeans_cm) vs log(n)
                    plt.plot(np.log(lambda_jeans_cm), x, lw=1, color='b')
                    plt.ylabel('log(n) [H/cc]')
                    plt.xlabel('log(lambda_Jeans) [cm]')
                    plt.grid(b=True, which='both', axis='both')
                    #plt.legend()
                    plt.title('log(lambda_Jeans) vs log(n)')
                    plt.savefig('log(lambda_Jeans)vsX_fromcode.png'.format(i=i))
                    plt.clf()
                    #plotting X_H2 vs log(n_H)
                    plt.plot(x, X_H2, lw=1, color='b')
                    plt.xlabel('log(n_H) [H/cc]')
                    plt.ylabel('X_H2')
                    plt.grid(b=True, which='both', axis='both')
                    #plt.legend()
                    plt.title('log(n_H) vs X_H2')
                    plt.savefig('log(n_H)vsX_H2.png'.format(i=i))
                    plt.clf()

if __name__=='__main__':
    #myobject = H2_density()
    H2_density.plotting()


'''
'''
#----- on 18.11.18 -----
# ---------------
# Code attempting to create multi-dimensional array for plotting
# ---------------

#n_H_mean_arr = np.array([1e-3, 1e-2, 1e-1, 1e1, 1e2, 1e3])
n_H_range = 100
x_mean = 1
pdf = np.zeros([6, 100])
x = np.zeros([6, 100])
s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
lambda_jeans_cm = np.zeros([6, 100])    # cm
tau = np.zeros([6, 100])     #optical depth
n_LW = np.zeros([6, 100])    #number of Lyman-Werner photons
n_H2 = np.zeros([6, 100])
X_H2 = np.zeros([6, 100])
n_H_bar = np.zeros([6, 100])
n_H2_bar = np.zeros([6, 100])
tot_n_H_bar = np.zeros([6, 100])
tot_n_H2_bar = np.zeros([6, 100])
X_H2_bar = np.zeros([6, 100])
tot_X_H2_bar = np.zeros([6, 100])

for n_H_mean_ctr in range(100):
    n_H_mean = (np.logspace(-4, 5, 100) * n_H_mean_ctr) # [H] cm^-3
    for i in range(1, n_H_range):
        n_H[n_H_mean_ctr][i] = [n_H_mean_ctr][(np.logspace(-4, 5, 100) * i)] # [H] cm^-3
        x[n_H_mean_ctr][i] = np.log(n_H[n_H_mean_ctr][i]/n_H_mean)
        pdf[n_H_mean_ctr][i] = make_pdf(x[n_H_mean_ctr][i], x_mean, s)
        lambda_jeans_cm[n_H_mean_ctr][i] = calc_lambda_jeans(T_mean, n_H[n_H_mean_ctr][i])
        tau[n_H_mean_ctr][i] = calc_optical_depth(n_H[n_H_mean_ctr][i], lambda_jeans_cm[n_H_mean_ctr][i])
        n_LW[n_H_mean_ctr][i] = calc_num_LW(tau[n_H_mean_ctr][i])
        X_H2[n_H_mean_ctr][i] = calc_X_H2(n_H[n_H_mean_ctr][i], n_LW[n_H_mean_ctr][i], Z)
        n_H2[n_H_mean_ctr][i] = calc_num_H2(n_H[n_H_mean_ctr][i], X_H2[n_H_mean_ctr][i])
        n_H2_bar[n_H_mean_ctr][i] = n_H[n_H_mean_ctr][i] * pdf[n_H_mean_ctr][i] * X_H2[n_H_mean_ctr][i]
        for k, item1 in enumerate(n_H[n_H_mean_ctr][i]):
            tot_n_H2_bar[n_H_mean_ctr][k] += quad(calc_n_H2_bar, -np.inf, np.inf, args=(n_H2_bar[n_H_mean_ctr][k], ))[0, 0]
        X_H2_bar[n_H_mean_ctr][i] =  n_H_mean/n_H2_bar[n_H_mean_ctr][i]
        tot_X_H2_bar[n_H_mean_ctr][i] = calc_X_H2_bar(n_H_mean, tot_n_H2_bar[n_H_mean_ctr][i])


#-----------------------

'''

'''
# ---------------
# Code for testing scipy.integrate.quad
# ---------------
def f(x):
    return x**2

if __name__=='__main__':
    f = quad(f, 0, 1)
    print (f)
'''
