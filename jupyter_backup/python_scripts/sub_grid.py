# This is code just for the sub_grid model, and
# to save X_H2, X_CO, n_CO, lambda_jeans arrays
# to use in rad_transfer code later.
# Values from simulation are not used here, this is purely theoretical.
# pynbody is only imported for assigning units to different quantities, although
# after importing the saved array the units are not imported for some reason, so 
# importing pynbody is useless. While imorting, pynbody can be used to assign the units again.


# Also, X_H2_bar and X_CO_bar are calculated and stored in:
#     reference_plots => for purely theoretical values
#     3.8.py, 4.0.py => for 5kpc and 15kpc versions by using simulation values
    

import numpy as np
import pynbody

def make_pdf(s, s_bar, sigma_s):
    pdf = (1./np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_bar)/sigma_s)**2)))
    return pdf

def calc_lambda_jeans(n_H, T_mean, m_p, K_b):
    lambda_jeans = (np.sqrt(K_b * T_mean/m_p) / np.sqrt(4* np.pi * G * n_H * m_p))
    return lambda_jeans

def calc_n_LW(n_H, G_o, lambda_jeans, Z, m_p):
    kappa = 1000 * m_p * Z
    rad_field_outside = G_o #in solar units
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    n_LW = rad_field_outside * exp_tau
    return n_LW

def calc_X_H2(n_H, Z, n_LW):
    DC = pynbody.array.SimArray(1.7e-11, "cm**2 g**-1 s**-1")
    CC = pynbody.array.SimArray(2.5e-17, "cm**3 s**-1")
    numerator = DC * n_LW
    denominator = CC * Z * n_H
    X_H2 = 1. / (2 + (numerator/denominator) )
    return X_H2

def calc_n_LW_ss(n_H, n_H2, G_o, lambda_jeans, m_p):
    kappa = 1000 * m_p * Z
    rad_field_outside = G_o #in solar units
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    N_H2 = n_H2 * lambda_jeans
    term1 = pynbody.array.SimArray((0.965/((1+(N_H2/5e14))**2)), "1")
    term2 = ( (0.035/np.sqrt(1+(N_H2/5e14))) * np.exp(-1*np.sqrt(1+(N_H2/5e14))/1180) )
    S_H2 = term1 + term2
    n_LW_ss = rad_field_outside * exp_tau * S_H2
    return n_LW_ss

def self_shielding_iterations(n_H, G_o, lambda_jeans, Z, m_p):
    n_LW = np.zeros(100)
    X_H2 = np.zeros(100)
    n_H2 = np.zeros(100)
    n_LW_ss = np.zeros(100)
    S_H2_ss = np.zeros(100)
    N_H2_ss = np.zeros(100)
    X_H2_ss = np.zeros(100)
    n_H2_ss = np.zeros(100)
    ctr = 16
    i = 0
    n_LW = calc_n_LW(n_H, G_o, lambda_jeans, Z, m_p)
    X_H2 = calc_X_H2(n_H, Z, n_LW)
    n_H2 = n_H * X_H2
    n_H2_ss = n_H2
    while i<ctr:
        n_LW_ss = calc_n_LW_ss(n_H, n_H2_ss, G_o, lambda_jeans, m_p)
        X_H2_ss = calc_X_H2(n_H, Z, n_LW_ss)
        n_H2_ss = n_H * X_H2_ss
        i += 1
    return n_LW, n_H2, X_H2, n_LW_ss, X_H2_ss, n_H2_ss

def calc_X_CO(n_H, n_H2, n_LW):
    rate_CHX = pynbody.array.SimArray(5.0e-10 * n_LW, "1")
    rate_CO = pynbody.array.SimArray(1.0e-10 * n_LW, "1") 
    x0 = pynbody.array.SimArray(2.0e-4, "1")
    k0 = pynbody.array.SimArray(5.0e-16, "1") #cm**3 s**-1
    k1 = pynbody.array.SimArray(5.0e-10, "1") #cm**3 s**-1
    factor_beta = pynbody.array.SimArray(rate_CHX/(n_H*k1*x0), "1")
    beta = 1./(1.+factor_beta)
    factor_CO = pynbody.array.SimArray(rate_CO/(n_H2*k0*beta), "1")
    X_CO = 1./(1.+factor_CO)
    return X_CO

def calc_n_CO(n_H, X_CO, Z):
    abundance_Ctot = 1e-4 # n_C/n_H as defined by nucleosynthesis
    return n_H * abundance_Ctot * Z * X_CO # CO/cc

def CO_density(mach_no, n_H_mean, Z, G_o, T_mean, m_p, K_b):
    s = np.zeros(100)
    pdf = np.zeros(100)
    lambda_jeans = np.zeros(100)
    n_H = np.zeros(100)
    X_CO = np.zeros(100)
    n_CO = np.zeros(100)

    sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -7*sigma_s + s_bar
    smax = 7*sigma_s + s_bar
    ds = (smax - smin)/100

    for i in range(0, 100):
        s[i] = smin + i*ds

    n_H = n_H_mean * np.exp(s)
    pdf = make_pdf(s, s_bar, sigma_s)
    lambda_jeans = calc_lambda_jeans(n_H, T_mean, m_p, K_b)
    
#     if (T>=1e4) | (n_H_mean <= 1e-2) :
#         n_LW = 0
#         X_H2 = 0
#         n_H2 = 0
#         n_LW_ss = 0
#         S_H2_ss = 0
#         N_H2_ss = 0
#         X_H2_ss = 0
#         n_H2_ss = 0
#         X_H2_bar = 0
#         X_CO = 0
#         n_CO = 0
#         X_CO_bar = 0
#     else:
    n_LW, n_H2, X_H2, n_LW_ss, X_H2_ss, n_H2_ss = self_shielding_iterations(n_H, G_o, lambda_jeans, Z, m_p)

    X_CO = calc_X_CO(n_H, n_H2, n_LW)
    n_CO = calc_n_CO(n_H, X_CO, Z)

    return n_H, n_H2, X_H2, n_LW, n_LW_ss, X_H2_ss, n_H2_ss, X_CO, n_CO, pdf, lambda_jeans

if __name__=='__main__':
    m_p = pynbody.array.SimArray(1.672621777e-24, "g")
    K_b = pynbody.array.SimArray(1.38064852e-16, "cm**2 g s**-2 K**-1")
    G = pynbody.array.SimArray(6.67259e-8, "cm**3 g**-1 s**-2")
    T_mean = pynbody.array.SimArray(10., "K")
    mach_no = pynbody.array.SimArray(10., "1")
    Z = 0.02/0.02
    G_o = 1
    n_H_mean = pynbody.array.SimArray(100., "cm**-3")

    n_H, n_H2, X_H2, n_LW, n_LW_ss, X_H2_ss, n_H2_ss, X_CO, n_CO, pdf, lambda_jeans = CO_density(mach_no, 
                                                                                                 n_H_mean, Z, G_o, 
                                                                                                 T_mean, m_p, K_b)
    
    np.save('outputs/sub_grid/n_H.npy', n_H)
    np.save('outputs/sub_grid/n_H2.npy', n_H2)
    np.save('outputs/sub_grid/X_H2.npy', 2*X_H2)
    np.save('outputs/sub_grid/n_LW.npy', n_LW)
    np.save('outputs/sub_grid/n_LW_ss.npy', n_LW_ss)
    np.save('outputs/sub_grid/X_H2_ss.npy', 2*X_H2_ss)
    np.save('outputs/sub_grid/n_H2_ss.npy', n_H2_ss)
    np.save('outputs/sub_grid/X_CO.npy', X_CO)
    np.save('outputs/sub_grid/n_CO.npy', n_CO)
    np.save('outputs/sub_grid/pdf.npy', pdf)
    np.save('outputs/sub_grid/lambda_jeans.npy', lambda_jeans)
    
    