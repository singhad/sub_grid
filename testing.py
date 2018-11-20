'''# ---------------
# Code for plotting X_H2 vs log(n_H) for different values of Z, mach_no, and
# n_H_mean
# ---------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import os

#function to generate the PDF
def make_pdf(x, x_mean, vel_disp):
    pdf = (1/np.sqrt(2*np.pi*(vel_disp**2))) * (np.exp(-0.5*(((x - x_mean)/vel_disp)**2)))
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

def calc_n_H2(n_H, X_H2):
    n_H2 = n_H * X_H2   #H2/cc
    return n_H2

def plot_X_H2_vs_x(X_H2, x, label, value):
    plt.plot(x, X_H2, label =(label + '= ' + str(value)))
    plt.xlabel('log(n_H)')
    plt.ylabel('X_H2')
    plt.legend()
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs X_H2 - varying ' + label)

def plotting(x, pdf, lambda_jeans_cm, n_H, X_H2, X_H2_bar, tot_n_H_bar, tot_n_H2_bar):
    path = 'multiple_plots'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    #plotting log(n) vs log(PDF)
    plt.plot(x, np.log(pdf))
    plt.xlabel('log(n) [H/cc]')
    plt.ylabel('log(PDF)')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n) vs log(PDF)')
    #plt.savefig(os.path.join('XvslogPDF','image_%04d.png'.format()))
    plt.savefig('XvslogPDF.png'.format())
    plt.clf()

    #plotting log(n) vs PDF
    plt.plot(x, pdf)
    plt.xlabel('log(n) [H/cc]')
    plt.ylabel('PDF')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n) vs PDF')
    #plt.savefig(os.path.join('XvsPDF','image_%04d.png'.format()))
    plt.savefig('XvsPDF.png'.format())
    plt.clf()

    #plotting log of Jeans length in cm log(lambda_jeans_cm) vs log(n)
    plt.plot(np.log(lambda_jeans_cm), x)
    plt.ylabel('log(n) [H/cc]')
    plt.xlabel('log(lambda_Jeans) [cm]')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(lambda_Jeans) vs log(n)')
    #plt.savefig(os.path.join('log(lambda_Jeans)vsX','image_%04d.png'.format()))
    plt.savefig('log(lambda_Jeans)vsX.png'.format())
    plt.clf()

    #plotting X_H2 vs log(n_H)
    plt.plot(x, X_H2)
    plt.xlabel('log(n_H) [H/cc]')
    plt.ylabel('X_H2')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs X_H2')
    #plt.savefig(os.path.join('log(n_H)vsX_H2','image_%04d.png'.format()))
    plt.savefig('log(n_H)vsX_H2.png'.format())
    plt.clf()

if __name__=='__main__':
    path = 'multiple_plots'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    m_p = 1.672621777e-24  # g
    T_mean = 10.    #K
    mach_no_arr = np.array([1., 5., 10., 15., 20., 50., 100.])
    Z_arr = np.array([0.001, 0.01, 0.1, 1., 10., 100.])
    n_H_mean_arr = np.array([1e-1, 1, 1e1, 1e2, 1e3, 1e4])
    n_H_range = 100
    x_mean = 1
    pdf = np.zeros(100)
    x = np.zeros(100)
    lambda_jeans_cm = np.zeros(100)    # cm
    tau = np.zeros(100)     #optical depth
    n_LW = np.zeros(100)    #number of Lyman-Werner photons
    n_H2 = np.zeros(100)
    X_H2 = np.zeros(100)

    choice = [1, 2, 3]
    for ch in choice:
        value = 0
        label = ""
        if ch == 1:
            mach_no = 5
            vel_disp = np.sqrt(np.log(1 + ((0.3 * mach_no)**2))) #vel_disp in pdf
            n_H_mean = 100
            label = "Metallicity "
            for i in range(0, len(Z_arr)):
                Z = Z_arr[i]
                value = Z
                for l in range(1, n_H_range):
                    n_H = (np.logspace(-4, 5, 100) * l) # [H] cm^-3
                    x = np.log(n_H/n_H_mean)
                    pdf = make_pdf(x, x_mean, vel_disp)
                    lambda_jeans_cm = calc_lambda_jeans(T_mean, n_H)
                    tau = calc_optical_depth(n_H, lambda_jeans_cm)
                    n_LW = calc_num_LW(tau)
                    X_H2 = calc_X_H2(n_H, n_LW, Z)
                    n_H2 = calc_n_H2(n_H, X_H2)
                plot_X_H2_vs_x(X_H2, x, label, value)
            plt.savefig('XvslogPDF-diff Z.png'.format())
            plt.clf()

        elif ch == 2:
            Z = 1
            n_H_mean = 100
            label = "Mach no "
            for i in range(0, len(mach_no_arr)):
                mach_no = mach_no_arr[i]
                vel_disp = np.sqrt(np.log(1 + ((0.3 * mach_no)**2))) #vel_disp in pdf
                value = mach_no
                for l in range(1, n_H_range):
                    n_H = (np.logspace(-4, 5, 100) * l) # [H] cm^-3
                    x = np.log(n_H/n_H_mean)
                    pdf = make_pdf(x, x_mean, vel_disp)
                    lambda_jeans_cm = calc_lambda_jeans(T_mean, n_H)
                    tau = calc_optical_depth(n_H, lambda_jeans_cm)
                    n_LW = calc_num_LW(tau)
                    X_H2 = calc_X_H2(n_H, n_LW, Z)
                    n_H2 = calc_n_H2(n_H, X_H2)
                plot_X_H2_vs_x(X_H2, x, label, value)
            plt.savefig('XvslogPDF-diff Mach no.png'.format())
            plt.clf()


        elif ch == 3:
            mach_no = 5
            vel_disp = np.sqrt(np.log(1 + ((0.3 * mach_no)**2))) #vel_disp in pdf
            Z = 1
            label = "n_H_mean "
            for i in range(0, len(n_H_mean_arr)):
                n_H_mean = n_H_mean_arr[i] # [H] cm^-3
                value = n_H_mean
                for l in range(1, n_H_range):
                    n_H = (np.logspace(-4, 5, 100) * l) # [H] cm^-3
                    x = np.log(n_H/n_H_mean)
                    pdf = make_pdf(x, x_mean, vel_disp)
                    lambda_jeans_cm = calc_lambda_jeans(T_mean, n_H)
                    tau = calc_optical_depth(n_H, lambda_jeans_cm)
                    n_LW = calc_num_LW(tau)
                    X_H2 = calc_X_H2(n_H, n_LW, Z)
                    n_H2 = calc_n_H2(n_H, X_H2)
                plot_X_H2_vs_x(X_H2, x, label, value)
            plt.savefig('XvslogPDF-diff n_H_mean.png'.format())
            plt.clf()'''

#plotting(x, pdf, lambda_jeans_cm, n_H, X_H2, X_H2_bar, tot_n_H_bar, tot_n_H2_bar)


# ---------------
# Code using scipy to integrate n_H_bar to find tot_n_H_bar
# ---------------

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad
import os

#function to generate the PDF
def make_pdf(x, x_mean, vel_disp):
    pdf = (1/np.sqrt(2*np.pi*(vel_disp**2))) * (np.exp(-0.5*(((x - x_mean)/vel_disp)**2)))
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

def calc_n_H2(n_H, X_H2):
    n_H2 = n_H * X_H2   #H2/cc
    return n_H2

def calc_tot_n_H_bar(n_H, n_H_bar):
    return n_H_bar

def calc_tot_n_H2_bar(n_H, n_H2_bar):
    return n_H2_bar

def plotting(x, pdf, lambda_jeans_cm, n_H, X_H2, X_H2_bar, tot_n_H_bar, tot_n_H2_bar):
    path = 'with_integration'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    #plotting log(n) vs log(PDF)
    plt.plot(x, np.log(pdf))
    plt.xlabel('log(n) [H/cc]')
    plt.ylabel('log(PDF)')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n) vs log(PDF)')
    #plt.savefig(os.path.join('XvslogPDF','image_%04d.png'.format()))
    plt.savefig('XvslogPDF.png'.format())
    plt.clf()

    #plotting log(n) vs PDF
    plt.plot(x, pdf)
    plt.xlabel('log(n) [H/cc]')
    plt.ylabel('PDF')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n) vs PDF')
    #plt.savefig(os.path.join('XvsPDF','image_%04d.png'.format()))
    plt.savefig('XvsPDF.png'.format())
    plt.clf()

    #plotting log of Jeans length in cm log(lambda_jeans_cm) vs log(n)
    plt.plot(np.log(lambda_jeans_cm), x)
    plt.ylabel('log(n) [H/cc]')
    plt.xlabel('log(lambda_Jeans) [cm]')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(lambda_Jeans) vs log(n)')
    #plt.savefig(os.path.join('log(lambda_Jeans)vsX','image_%04d.png'.format()))
    plt.savefig('log(lambda_Jeans)vsX.png'.format())
    plt.clf()

    #plotting X_H2 vs log(n_H)
    plt.plot(x, X_H2)
    plt.xlabel('log(n_H) [H/cc]')
    plt.ylabel('X_H2')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs X_H2')
    #plt.savefig(os.path.join('log(n_H)vsX_H2','image_%04d.png'.format()))
    plt.savefig('log(n_H)vsX_H2.png'.format())
    plt.clf()

    #plotting log(tot_n_H_bar) vs X_H2_bar
    plt.plot(np.log(tot_n_H_bar), X_H2_bar)
    plt.xlabel('log(tot_n_H_bar) [H/cc]')
    plt.ylabel('X_H2_bar')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(tot_n_H_bar) vs X_H2_bar')
    #plt.savefig(os.path.join('n_HvsX_H2_bar','image_%04d.png'.format()))
    plt.savefig(os.path.join('log(tot_n_H_bar)vsX_H2_bar.png'.format()))
    plt.clf()


if __name__=='__main__':
    m_p = 1.672621777e-24  # g
    T_mean = 10.    #K
    mach_no = 5
    Z = 1
    n_H_mean = 100
    n_H_range = 100
    x_mean = 1
    vel_disp = np.sqrt(np.log(1 + ((0.3 * mach_no)**2))) #vel_disp in pdf
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
    tot_n_H_bar = np.zeros(100)
    tot_n_H2_bar = np.zeros(100)

    for i in range(1, n_H_range):
        n_H = (np.logspace(-4, 5, 100) * i) # [H] cm^-3
        x = np.log(n_H/n_H_mean)
        pdf = make_pdf(x, x_mean, vel_disp)
        lambda_jeans_cm = calc_lambda_jeans(T_mean, n_H)
        tau = calc_optical_depth(n_H, lambda_jeans_cm)
        n_LW = calc_num_LW(tau)
        X_H2 = calc_X_H2(n_H, n_LW, Z)
        n_H2 = calc_n_H2(n_H, X_H2)
        n_H_bar = n_H * pdf
        n_H2_bar = n_H * pdf * X_H2
        for k, item in enumerate(n_H2_bar):
            tot_n_H2_bar[k] = quad(calc_tot_n_H2_bar, -np.inf, np.inf, args=(n_H2_bar[k], ))[0]
        for l, item in enumerate(n_H_bar):
            tot_n_H_bar[l] = quad(calc_tot_n_H_bar, -np.inf, np.inf, args=(n_H_bar[l], ))[0]
        X_H2_bar =  tot_n_H_bar/tot_n_H2_bar

    plotting(x, pdf, lambda_jeans_cm, n_H, X_H2, X_H2_bar, tot_n_H_bar, tot_n_H2_bar)
