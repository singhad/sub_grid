#!/usr/bin/env python

import timing
import numpy as np
import pynbody
from multiprocessing import Pool
from michaels_functions import (center_and_r_vir, remove_bulk_velocity,
                                read_unit_from_info)

def get_filename(species):
    # filename is already given
    if (species[-4:] == '.dat') or (species[-4:] == '.txt'):
        return species
    # molecule is chosen
    THIS_FOLDER = os.path.dirname(os.path.abspath(__file__))
    database = os.path.join(THIS_FOLDER, 'LAMDA')
    if species == 'HCO+':
        filename = os.path.join(database, 'HCO+.dat')
    elif species == 'H13CO+':
        filename = os.path.join(database, 'H13CO+.dat')
    elif species == 'N2H+':
        filename = os.path.join(database, 'N2H+.dat')
    elif species == 'SiO':
        filename = os.path.join(database, 'SiO.dat')
    elif species == 'HNC':
        filename = os.path.join(database, 'HNC.dat')
    elif species == 'HCN':
        filename = os.path.join(database, 'HCN.dat')
    elif species == 'CO':
        filename = os.path.join(database, 'CO.dat')
    else:
        print('Unknow species. Chose from HCO+, H13CO+, N2H+, SiO, HNC, HCN, CO')
        print('or provide a LAMDA datafile.')
        exit()

    return filename


def read_file(species, c_cgs, h_ev, K_b_ev):
    filename = get_filename(species)
    f = open(filename, 'r')

    f.readline()
    species = f.readline()

    f.readline()
    mu = float(f.readline())  # molecular weight

    f.readline()
    num_lvls = int(f.readline())  # number of energy levels

    # read energy levels: energy E, statistical weight g
    f.readline()
    E = []
    g = []
    for l in range(num_lvls):
        words = f.readline().split()
        E.append(float(words[1]) *c_cgs*h_ev) # cm^-1 -> eV
        g.append(float(words[2]))

    f.readline()
    num_trans = int(f.readline())  # number of radiative transistions

    # read transistions: upper lvl, lower lvl, A-coefficient, frequency
    f.readline()
    A = np.zeros((num_lvls, num_lvls))
    freq = np.zeros((num_lvls, num_lvls))
    for t in range(num_trans):
        words = f.readline().split()
        up = int(words[1]) - 1
        low = int(words[2]) - 1
        if up-low==1:
            A[up][low] = float(words[3])  # s^-1
            freq[up][low] = float(words[4]) * 1e9  # GHz -> Hz
            #freq[low][up] = freq[up][low] #un-comment this only if low->up transitions are also allowed

    # compute B-coefficient via Einstein relations
    # Bij = coeff for stimulated emission, Bji = coeff for extinction (j<i)
    B = np.zeros((num_lvls, num_lvls))
    for i in range(0, num_lvls):
        for j in range(0, i):
            if A[i][j] != 0:
                B[i][j] = A[i][j] * (c_cgs**2) / (2*h_ev * (freq[i][j])**3) # cm2/(eV*s)
                B[j][i] = B[i][j] * g[i]/g[j]

    f.close()
    E = np.array(E) #eV
    g = np.array(g) 
    freq = np.array(freq) #Hz
    A = np.array(A) #s-1
    B = np.array(B) #cm2/(eV*s)
    return mu, num_lvls, E, g, freq, A, B


def partion_function(T, num_lvls, g, E):
    Z=0.0
    for i in range(0,num_lvls):
        Z = Z + g[i]*np.exp(-E[i]/(K_b_ev*T))
    return np.asarray(Z)


def calc_lvlpops_partion(T, num_lvls, g, E):
    ni = np.zeros(num_lvls)
    Z = partion_function(T, num_lvls, g, E)
    for i in range(0, num_lvls):
        ni[i] = g[i]*np.exp(-E[i]/(K_b_ev*T)) / Z
    return np.asarray(ni), Z


def make_pdf(s, s_bar, sigma_s):
    return ((1./np.sqrt(2*np.pi*(sigma_s**2))) *
            (np.exp(-0.5*(((s - s_bar)/sigma_s)**2))))


def calc_lambda_jeans(n_H, T_mean, m_p, K_b):
    return np.asarray(np.sqrt(K_b * T_mean/m_p) /
                      np.sqrt(4.*np.pi * G * n_H * m_p))


def calc_n_LW(n_H, G_o, lambda_jeans, metal, m_p):
    kappa = 1000 * m_p * metal
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    return np.asarray(G_o * exp_tau)


def calc_X_H2(n_H, metal, n_LW):
    DC = 1.7e-11
    CC = 2.5e-17            # cm3 s-1
    numerator = DC * n_LW
    denominator = CC * metal * n_H
    X_H2 = 1. / (2. + (numerator/denominator))
    return np.asarray(X_H2)


def calc_n_LW_ss(n_H, n_H2, metal, G_o, lambda_jeans, m_p):
    kappa = 1000 * m_p * metal
    rad_field_outside = G_o  # in solar units
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    N_H2 = n_H2 * lambda_jeans
    term1 = 0.965/((1+(N_H2/5e14))**2)
    term2 = ((0.035/np.sqrt(1+(N_H2/5e14))) *
             np.exp(-1*np.sqrt(1+(N_H2/5e14))/1180))
    S_H2 = term1 + term2
    n_LW_ss = rad_field_outside * exp_tau * S_H2
    return np.asarray(n_LW_ss)


def self_shielding_iterations(n_H, G_o, lambda_jeans, metal, m_p):
    ctr = 16
    n_LW = calc_n_LW(n_H, G_o, lambda_jeans, metal, m_p)
    X_H2 = calc_X_H2(n_H, metal, n_LW)
    n_H2 = n_H * X_H2
    n_H2_ss = n_H2
    for _ in range(ctr):
        n_LW_ss = calc_n_LW_ss(n_H, n_H2_ss, metal, G_o, lambda_jeans, m_p)
        X_H2_ss = calc_X_H2(n_H, metal, n_LW_ss)
        """
        if (np.sum(np.square(n_H2_ss - n_H * X_H2_ss)) < 1e-5):
            n_H2_ss = n_H * X_H2_ss
            break
        """
        n_H2_ss = n_H * X_H2_ss
    return n_LW, n_H2, n_LW_ss, 2*X_H2_ss, 2*n_H2_ss


def calc_integral(s, pdf, X, ds):
    return np.sum(np.exp(s) * pdf * X * ds)


def calc_X_CO(n_H, n_H2, n_LW):
    rate_CHX = 5.0e-10 * n_LW
    rate_CO = 1.0e-10 * n_LW
    x0 = 2.0e-4
    k0 = 5.0e-16  # cm3 s-1
    k1 = 5.0e-10  # cm3 s-1
    factor_beta = rate_CHX/(n_H*k1*x0)
    beta = 1./(1.+factor_beta)
    factor_CO = rate_CO/(n_H2*k0*beta)
    X_CO = 1./(1.+factor_CO)
    return np.asarray(X_CO)


def calc_n_CO(n_H, X_CO, metal):
    abundance_Ctot = 1e-4 # n_C/n_H as defined by nucleosynthesis
    return n_H * abundance_Ctot * metal * X_CO # CO/cc


def tau_LVG(N, nu, lambda_jeans, x_1, x_0, B_10, B_01, c_s_CO, c_cgs, h_ev):
    # units: eV*s * Hz * cm * 1/cm3 * cm2/(eV*s) * 1/(Hz) = none
    delta_nu = c_s_CO*nu/c_cgs #"Hz"
    return h_ev*nu*lambda_jeans*N*((x_0*B_01)-(x_1*B_10)) / (4*np.pi*delta_nu)
 
    
def beta_LVG(tau):
    beta_lvg = np.zeros(len(tau))
    for i in range(0, len(tau)):
        if tau[i] < 0.01:
            beta_lvg[i] = 1. - tau[i]/2.
        elif tau[i] > 100.:
            beta_lvg[i] = 1./tau[i]
        else:
            beta_lvg[i] = (1.0 - np.exp(-tau[i])) / tau[i]
    return beta_lvg


def calc_integrated_emissivity(N, nu, x_1, A_10, h_ev):
    #units: eV*s * Hz * cm-3 * s-1 = eV/cm3/s
    j_10 = h_ev * nu * N * x_1 * A_10
    return j_10


def inside_loop(mach_no, n_H_mean, metal, G_o, T_mean, m_p, K_b, K_b_ev, 
                K_b_erg, c_cgs, h_ev, c_s_CO, E, A_10, B_10, B_01, nu, 
                x_1, x_0, eV_to_ergs, cell_width):
    sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -7*sigma_s + s_bar
    smax = 7*sigma_s + s_bar

    s = np.linspace(smin, smax, num=100, endpoint=False)
    ds = np.diff(s)[0]

    n_H = n_H_mean * np.exp(s)
    pdf = make_pdf(s, s_bar, sigma_s)
    lambda_jeans = calc_lambda_jeans(n_H, T_mean, m_p, K_b)

    n_LW, n_H2, n_LW_ss, X_H2_ss, n_H2_ss = self_shielding_iterations(
        n_H, G_o, lambda_jeans, metal, m_p)

    X_H2_bar = calc_integral(s, pdf, X_H2_ss, ds)
    X_CO = calc_X_CO(n_H, n_H2, n_LW)
    n_CO = calc_n_CO(n_H, X_CO, metal)
    X_CO_bar = calc_integral(s, pdf, X_CO, ds)
    
    tau_nu = tau_LVG(n_CO, nu, lambda_jeans, x_1, x_0, B_10, B_01, c_s_CO, c_cgs, h_ev)
    beta_nu = beta_LVG(tau_nu)
    j_10 = calc_integrated_emissivity(n_CO, nu, x_1, A_10, h_ev)  # eV/cm3/s
    m_H2 = m_p * (cell_width**3) * X_H2_ss  #g cm3
    l_CO = j_10 * beta_nu * eV_to_ergs * m_H2/m_p  # ergs/s
    
    alpha_CO = ((c_cgs/nu)**3)*1e-5*((3.24078e-19)**2)/(2*K_b_erg) #conversion factor: erg/s to K km s-1 pc2
    l_CO_SI = l_CO * alpha_CO # K km s-1 pc2
    
    m_H2_bar = np.sum(m_H2*n_H*pdf*ds)
    l_CO_bar = np.sum(l_CO*pdf*ds)
    l_CO_SI_bar = np.sum(l_CO_SI*pdf*ds)

    return X_H2_bar, X_CO_bar, l_CO_bar, l_CO_SI_bar, m_H2_bar


if __name__ == '__main__':
    run = "hydro_59"
    out = "output_00050"
    path = "bulk1/data_2/" + run + "/output/"
    data = pynbody.load(path + out)
    aexp = data.properties['a']
    data.physical_units()

    r_vir = center_and_r_vir(data, aexp, path)
    remove_bulk_velocity(data)
    r_e = 0.1 * r_vir

    sph_5 = pynbody.filt.Sphere(radius='%f kpc' % r_e)
    region = data[sph_5]

    omega_b, unit_l, unit_d, unit_t = read_unit_from_info(data)

    m_p = pynbody.array.SimArray(1.672621777e-24, "g")
    K_b = pynbody.array.SimArray(1.38064852e-16, "cm**2 g s**-2 K**-1")
    G = pynbody.array.SimArray(6.67259e-8, "cm**3 g**-1 s**-2")
    T_mean = pynbody.array.SimArray(10., "K")
    K_b_ev = pynbody.array.SimArray(8.617e-5, "eV K**-1")
    K_b_erg = pynbody.array.SimArray(1.38064852e-16, "erg K**-1")
    c_cgs = pynbody.array.SimArray(2.99792458e10, "cm s**-1")
    h_ev = pynbody.array.SimArray(4.13566770e-15, "eV s")
    h_erg = pynbody.array.SimArray(6.626196e-27, "erg s")
    mH_cgs = pynbody.array.SimArray(1.6737236e-24, 'g')  # hydrogen mass
    eV_to_ergs = pynbody.array.SimArray(1.60218e-12, "erg eV**-1")
    cell_width_arr = pynbody.array.SimArray(region.gas["smooth"]*3.086e+21, "cm")

    rho = region.gas["rho"].in_units("m_p cm**-3")
    mass = region.gas["mass"].in_units('Msol')
    turb = np.sqrt(region.g["turb"] * 2./3.) * unit_l / unit_t / 1e5
    turb = pynbody.array.SimArray(turb*1e5, units="cm s**-1")

    temperature = region.g["temp"]
    c_s_arr = np.sqrt(K_b * temperature / m_p)

    mach_no_sim = np.array(turb / c_s_arr)

    m_p_1 = pynbody.array.SimArray(1.0, pynbody.units.m_p)
    n_H_mean_sim = np.array(rho / m_p_1)

    metal_arr = region.g["metal"]/0.02
    G_o = 1
    
    mu, num_lvls, E, g, freq, A, B = read_file('CO.txt',c_cgs, h_ev, K_b_ev)
    ni, Z = calc_lvlpops_partion(T_mean, num_lvls, g, E)
    c_s_CO = np.sqrt(K_b * T_mean/(mH_cgs*mu))
    u = 1 #upper level
    l = 0 #lower level
    nu = freq[u][l]
    x_1 = ni[u]
    x_0 = ni[l]
    A_10 = A[u][l]
    B_10 = B[u][l]
    B_01 = B[l][u]    
    
    mask_relevant = np.logical_and(temperature < 1e4, n_H_mean_sim > 1e-2)
    original_order = np.arange(len(rho))
    rel_cells = original_order[mask_relevant]
    X_H2_bar_cells = np.zeros(len(rho))
    X_CO_bar_cells = np.zeros(len(rho))
    l_CO_bar_cells = np.zeros(len(rho))
    l_CO_SI_bar_cells = np.zeros(len(rho))
    m_H2_bar_cells = np.zeros(len(rho))

    def map_to_loop(cell):
        mach_no = mach_no_sim[cell]
        n_H_mean = n_H_mean_sim[cell]
        metal = metal_arr[cell]
        cell_width = cell_width_arr[cell]

        X_H2_bar, X_CO_bar, l_CO_bar, l_CO_SI_bar, m_H2_bar = inside_loop(mach_no, n_H_mean, metal, 
        G_o, T_mean, m_p, K_b, K_b_ev, K_b_erg, c_cgs, h_ev, c_s_CO, E, A_10, B_10, B_01, nu, x_1, x_0, eV_to_ergs, cell_width)
        return X_H2_bar, X_CO_bar, l_CO_bar, l_CO_SI_bar, m_H2_bar

    pool = Pool(16)
    X_list = pool.map(map_to_loop, rel_cells)
    pool.close()
    pool.join()

    for j, cell in enumerate(rel_cells):
        X_H2_bar, X_CO_bar, l_CO_bar, l_CO_SI_bar, m_H2_bar = X_list[j]
        X_H2_bar_cells[cell] = X_H2_bar
        X_CO_bar_cells[cell] = X_CO_bar
        l_CO_bar_cells[cell] = l_CO_bar
        l_CO_SI_bar_cells[cell] = l_CO_SI_bar
        m_H2_bar_cells[cell] = m_H2_bar

    np.save('outputs/parallel_code_pdf_RT/X_H2_bar_' + run + '_' + out + '.npy', X_H2_bar_cells)
    np.save('outputs/parallel_code_pdf_RT/X_CO_bar_' + run + '_' + out + '.npy', X_CO_bar_cells)
    np.save('outputs/parallel_code_pdf_RT/mach_no_arr_' + run + '_' + out + '.npy', mach_no_sim)
    np.save('outputs/parallel_code_pdf_RT/n_H_mean_arr_' + run + '_' + out + '.npy', n_H_mean_sim)
    np.save('outputs/parallel_code_pdf_RT/metal_arr_' + run + '_' + out + '.npy', metal_arr)
    np.save('outputs/parallel_code_pdf_RT/T_' + run + '_' + out + '.npy', temperature)
    np.save('outputs/parallel_code_pdf_RT/l_CO_bar_' + run + '_' + out + '.npy', l_CO_bar_cells)
    np.save('outputs/parallel_code_pdf_RT/l_CO_SI_bar_' + run + '_' + out + '.npy', l_CO_SI_bar_cells)
    np.save('outputs/parallel_code_pdf_RT/m_H2_bar_' + run + '_' + out + '.npy', m_H2_bar_cells)
    np.save('outputs/parallel_code_pdf_RT/mass_' + run + '_' + out + '.npy', mass)