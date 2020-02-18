# Here, Z is converted into solar units:
#     Z_sim is in mass fraction units
#     Z_arr = Z_sim/0.02
#     kappa = 1000*m_p*Z_arr  --> in n_LW & n_LW_ss
#     denominator = CC * Z * n_H  --> in X_H2
#     n_CO = 1e-4*n_H*X_CO*Z_arr

import timing
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import pynbody
from michaels_functions import center_and_r_vir, remove_bulk_velocity
from matplotlib.colors import LogNorm
from matplotlib.pyplot import figure

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
    return np.array(Z)

def calc_lvlpops_partion(T, num_lvls, g, E):
    ni = np.zeros(num_lvls)
    Z = partion_function(T, num_lvls, g, E)
    for i in range(0, num_lvls):
        ni[i] = g[i]*np.exp(-E[i]/(K_b_ev*T)) / Z
    return np.array(ni), Z

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
    DC = 1.7e-11
    CC = 2.5e-17            #cm3 s-1
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
    return n_LW, n_H2, n_LW_ss, X_H2_ss, n_H2_ss

def calc_X_CO(n_H, n_H2, n_LW):
    rate_CHX = 5.0e-10 * n_LW
    rate_CO = 1.0e-10 * n_LW
    x0 = 2.0e-4
    k0 = 5.0e-16 #cm3 s-1
    k1 = 5.0e-10 #cm3 s-1
    factor_beta = rate_CHX/(n_H*k1*x0)
    beta = 1./(1.+factor_beta)
    factor_CO = rate_CO/(n_H2*k0*beta)
    X_CO = 1./(1.+factor_CO)
    return X_CO

def calc_n_CO(n_H, X_CO, Z):
    abundance_Ctot = 1e-4 # n_C/n_H as defined by nucleosynthesis
    return n_H * abundance_Ctot * Z * X_CO # CO/cc

def calc_line_profile(nu, c_s_CO, c_cgs):
    delta_nu = c_s_CO*nu/c_cgs #"Hz"
    return 1/delta_nu #Hz-1

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

def B_nu_ev(nu, temp, h_ev, K_b_ev, c_cgs):
    if nu==0.:
        return 0.
    if nu>0:
        x = h_ev*nu/(K_b_ev*temp) #units: none
        #units: eV*s * Hz3 / (cm2/s2) = eV * s3 * Hz3 * cm-2 = eV/s/Hz/cm2
        B_nu_ev = 2.0*h_ev*(nu**3) / ((c_cgs**2) * (np.exp(x)-1.0))
        return B_nu_ev

def calc_integrated_emissivity(N, nu, x_1, A_10, h_ev):
    #units: eV*s * Hz * cm-3 * s-1 = eV/cm3/s
    j_10 = h_ev * nu * N * x_1 * A_10
    return j_10

def inside_loop(mach_no, n_H_mean, metal, G_o, T_mean, T, G, m_p, K_b, K_b_ev, K_b_erg, c_cgs, h_ev, T_bg, eV, s, n_H, pdf,
                lambda_jeans, X_CO, n_CO, sigma_s, s_bar, smin, smax, ds, c_s_CO, freq, E, A_10, B_10, B_01, ni, Z, tau_nu,
                beta_nu, B_nu, j_10, phi_ij, nu, x_1, x_0, eV_to_ergs, cell_width):
    
    sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -7*sigma_s + s_bar
    smax = 7*sigma_s + s_bar
    ds = (smax - smin)/100
    l_CO_bar = 0.0
    l_CO_SI_bar = 0.0
    for i in range(0, 100):
        s[i] = smin + i*ds
    n_H = n_H_mean * np.exp(s)
    pdf = make_pdf(s, s_bar, sigma_s)
    lambda_jeans = calc_lambda_jeans(n_H, T_mean, m_p, K_b)
    
    if (T>=1e4) | (n_H_mean <= 1e-2) :
        n_LW = 0
        X_H2 = 0
        n_H2 = 0
        n_LW_ss = 0
        S_H2_ss = 0
        N_H2_ss = 0
        X_H2_ss = 0
        n_H2_ss = 0
        X_H2_bar = 0
        X_CO = 0
        n_CO = 0
        X_CO_bar = 0
    else:
        n_LW, n_H2, n_LW_ss, X_H2_ss, n_H2_ss = self_shielding_iterations(n_H, G_o, lambda_jeans, metal, m_p)
        X_CO = calc_X_CO(n_H, n_H2, n_LW)
        n_CO = calc_n_CO(n_H, X_CO, metal)
        X_H2_bar = 2 * np.sum(s*pdf*X_H2_ss*ds)
        X_CO_bar = np.sum(s*pdf*X_CO*ds)
        
    #B_nu = B_nu_ev(nu, temp, h_ev, K_b_ev, c_cgs)
    N = n_CO
    #Ni = ni*N
    tau_nu = tau_LVG(N, nu, lambda_jeans, x_1, x_0, B_10, B_01, c_s_CO, c_cgs, h_ev)
    beta_nu = beta_LVG(tau_nu)
    j_10 = calc_integrated_emissivity(N, nu, x_1, A_10, h_ev)  # eV/cm3/s
    l_CO = j_10 * beta_nu * eV_to_ergs  # ergs/cm3/s
    alpha_CO = ((c_cgs/nu)**3)*1e-5*((3.24078e-19)**2)/ (2*K_b_erg) #conversion factor: eV/cm3/s to K km s-1 pc2
    l_CO_SI = l_CO * alpha_CO # K km s-1 pc2
    
    l_CO_bar = np.sum(l_CO*pdf*ds)
    l_CO_SI_bar = np.sum(l_CO_SI*pdf*ds)
        
    return X_H2_bar, X_CO_bar, l_CO_bar, l_CO_SI_bar


if __name__=='__main__':
    path = "bulk1/data_2/hydro_59/output/"
    data = pynbody.load(path + "output_00050")
    aexp = data.properties['a']
    data.physical_units()
    r_vir = center_and_r_vir(data, aexp, path)
    remove_bulk_velocity(data)
    r_e = 0.1 * r_vir
    sph_5 = pynbody.filt.Sphere(radius = '%f kpc' %(r_e*1.4))
    region = data[sph_5]
    rho = region.gas["rho"].in_units("m_p cm**-3")
    f = open(data.filename + "/info_"+data.filename[-5:]+".txt","r")
    lines = f.readlines()
    f.close()
    for line in lines:
        if line[0:13]=="unit_l      =":
            print line[:-1]
            unit_l = float(line[14:-1])
        if line[0:13]=="unit_d      =":
            print line[:-1]
            unit_d = float(line[14:-1])
        if line[0:13]=="unit_t      =":
            print line[:-1]
            unit_t = float(line[14:-1])
        if line[0:13]=="omega_b     =":
            print line[:-1]
            omega_b = float(line[14:-1])
    
    m_p = pynbody.array.SimArray(1.672621777e-24, "g")
    K_b = pynbody.array.SimArray(1.38064852e-16, "cm**2 g s**-2 K**-1")
    G = pynbody.array.SimArray(6.67259e-8, "cm**3 g**-1 s**-2")
    T_mean = pynbody.array.SimArray(10., "K")
    K_b_ev = pynbody.array.SimArray(8.617e-5, "eV K**-1")
    K_b_erg = pynbody.array.SimArray(1.38064852e-16, "erg K**-1")
    c_si = pynbody.array.SimArray(2.99792458e5, "km s**-1")
    c_cgs = pynbody.array.SimArray(2.99792458e10, "cm s**-1")
    h_ev = pynbody.array.SimArray(4.13566770e-15, "eV s")
    h_erg = pynbody.array.SimArray(6.626196e-27, "erg s")
    mH_cgs = pynbody.array.SimArray(1.6737236e-24, 'g')  # hydrogen mass
    T_bg = pynbody.array.SimArray(2.73, "K")
    eV = pynbody.array.SimArray(6.241509e18, "J")
    L_sun = pynbody.array.SimArray(2.43418864387146974e+45, "eV s**-1")
    cell_width = pynbody.array.SimArray(region.gas["smooth"]*3.086e+21, "cm")
    M_sun = pynbody.array.SimArray(2e33, "g")
    eV_to_ergs = pynbody.array.SimArray(1.60218e-12, "erg eV**-1")
    cell_width_arr = pynbody.array.SimArray(region.gas["smooth"]*3.086e+21, "cm")
    
    turb = np.sqrt( region.g["turb"] * 2./3. ) * unit_l / unit_t / 1e5
    turb = pynbody.array.SimArray(turb*1e5, units = "cm s**-1")
    #turb_SI = pynbody.array.SimArray(turb, units = "km s**-1")

    temperature = region.g["temp"]
    c_s_arr = np.sqrt(K_b * temperature / m_p)

    mach_no_sim = turb / c_s_arr
    region.g["mach"] = mach_no_sim.in_units("1")

    m_p_1 = pynbody.array.SimArray(1.0, pynbody.units.m_p)
    n_H_mean_sim = rho / m_p_1

    metal_arr = region.g["metal"]/0.02
    G_o = 1

    mach_no_arr = mach_no_sim
    n_H_mean_arr = n_H_mean_sim
    mu, num_lvls, E, g, freq, A, B = read_file('CO.txt',c_cgs, h_ev, K_b_ev)
    ni, Z = calc_lvlpops_partion(T_mean, num_lvls, g, E)
    c_s_CO = np.sqrt(K_b * T_mean/(mH_cgs*mu))
    
    s = np.zeros(100)
    n_H = np.zeros(100)
    pdf = np.zeros(100)
    lambda_jeans = np.zeros(100)
    X_CO = np.zeros(100)
    n_CO = np.zeros(100)
    integral1 = 0.0
    integral2 = 0.0
    sigma_s = 0.0
    s_bar = 0.0
    smin = 0.0
    smax = 0.0
    ds = 0.0
    tau_nu = np.zeros(100)
    beta_nu = np.zeros(100)
    B_nu = np.zeros(100)
    j_10 = np.zeros(100)
    u = 1 #upper level
    l = 0 #lower level
    nu = freq[u][l]
    x_1 = ni[u]
    x_0 = ni[l]
    A_10 = A[u][l]
    B_10 = B[u][l]
    B_01 = B[l][u]
    phi_ij = calc_line_profile(nu, c_s_CO, c_cgs)

    X_H2_bar = np.zeros(len(n_H_mean_arr))
    X_CO_bar = np.zeros(len(n_H_mean_arr))
    l_CO_bar = np.zeros(len(n_H_mean_arr))
    l_CO_SI_bar = np.zeros(len(n_H_mean_arr))

    for m in range(0, len(mach_no_arr)):
        mach_no = mach_no_arr[m]
        n_H_mean = n_H_mean_arr[m]
        metal = metal_arr[m]
        T = temperature[m]
        cell_width = cell_width_arr[m]
        X_H2_bar[m], X_CO_bar[m], l_CO_bar[m], l_CO_SI_bar[m] = inside_loop(mach_no, n_H_mean, metal, G_o, T_mean, T, G, m_p, 
                                                            K_b, K_b_ev, K_b_erg, c_cgs, h_ev, T_bg, eV, s, n_H, pdf,
                                                            lambda_jeans, X_CO, n_CO, sigma_s, s_bar, smin, smax, ds, c_s_CO,
                                                            freq, E, A_10, B_10, B_01, ni, Z, tau_nu, beta_nu, B_nu, j_10,                                                                   phi_ij, nu, x_1, x_0, eV_to_ergs, cell_width)

    np.save('outputs/latest_sim_run/X_H2_bar.npy', X_H2_bar)
    np.save('outputs/latest_sim_run/X_CO_bar.npy', X_CO_bar)
    np.save('outputs/latest_sim_run/mach_no_arr.npy', mach_no_arr)
    np.save('outputs/latest_sim_run/n_H_mean_arr.npy', n_H_mean_arr)
    np.save('outputs/latest_sim_run/metal_arr.npy', metal_arr)
    np.save('outputs/latest_sim_run/T.npy', temperature)
    np.save('outputs/latest_sim_run/l_CO_bar.npy', l_CO_bar)
    np.save('outputs/latest_sim_run/l_CO_SI_bar.npy', l_CO_SI_bar)