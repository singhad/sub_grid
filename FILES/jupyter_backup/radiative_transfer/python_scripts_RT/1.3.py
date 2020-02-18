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

    # number of collision partners in the data file
    f.readline()
    num_partners = int(f.readline())

    C_all = []
    temps_all = []
    for partner in range(num_partners):
        # reference
        f.readline()
        line = f.readline()

        # number of collisional transitions
        f.readline()
        num_coll = int(f.readline())

        # number of temperatures in the table
        f.readline()
        num_temps = int(f.readline())

        # read the temperature values
        f.readline()
        words = f.readline().split()
        temps = np.zeros(num_temps)
        for t in range(num_temps):
            temps[t] = float(words[t])
            temps_all.append(temps)  # K

        # read collision coeff data: upper lvl, lower lvl, C-coefficient for each temp
        C = np.zeros((num_temps, num_lvls, num_lvls))
        f.readline()
        for col in range(num_coll):
            words = f.readline().split()
            i = int(words[1]) - 1
            j = int(words[2]) - 1
            for t in range(num_temps):
                C[t][i][j] = float(words[3+t]) #* 1.e-6 # cm3/s -> m3/s

        # calculate the inverse coefficient via LTE relation
        for i in range(num_lvls):
            for j in range(i):
                for t in range(num_temps):
                    if C[t][i][j] != 0:
                        C[t][j][i] = C[t][i][j] * np.exp(-(E[i]-E[j])/(K_b_ev*temps[t]))*g[i]/g[j]

        # add collision partner data to global array
        C_all.append(C)

    f.close()
    C_all = np.array(C_all) #cm3/s
    temps_all = np.array(temps_all) #K
    E = np.array(E) #eV
    g = np.array(g) 
    freq = np.array(freq) #Hz
    A = np.array(A) #s-1
    B = np.array(B) #cm2/(eV*s)
    return mu, num_lvls, E, g, freq, A, B, C_all, num_partners, temps_all, num_temps, num_coll

''' Load preset abundances and fraction for collision coefficients
    PARAMS:
      species = string with the particle name
    RETURNS:
      comp_fracs = list with the fraction of total collision partner density for each partner
      abundance = overall abundance of the molecule (assume n_mol = abundance*rho everywhere)'''
def load_species_info(species):

    if species == 'HCO+':
        comp_fracs = [1.0] # only 1 collision partner in HCO data
        abundance = 1.e-09 # = N_species/N_h2
    elif species == 'H13CO+':
        comp_fracs = [1.0] # only 1 collision partner in HCO data
        abundance = 2.e-11
    elif species == 'N2H+':
        comp_fracs = [1.0] # only 1 collision partner in HCO data
        abundance = 1.e-10
    elif species == 'SiO': # is seems to be unusually slow
        comp_fracs = [1.0] # only 1 collision partner in HCO data
        abundance = 7.7e-12
    elif species == 'HNC':
        comp_fracs = [1.0] # only 1 collision partner in HCO data
        abundance = 3.1e-10
    elif species == 'HCN':
        comp_fracs = [1.0, 0.0] # H2 and e
        abundance = 3.1e-11
    elif species == 'CO':
        comp_fracs = [0.66, 0.33] # para H2 and orhto H2
        #abundance = 1.e-4
        abundance = 1.0 # dummy for filenames
    else:
        print 'ERROR: Unsupported species'
        exit()

    return comp_fracs, abundance

''' Calculate net collision coeff for a gas consisting of different components at temp T
    Interpolates table betweem T values
    PARAMS:
      T = temperature (K)
      comp_fractions = fraction of the total density in each component
    RETRUN:
      C = netto collision coeff C[i][j] (m3/s) '''
def calc_total_C(T, comp_fractions):
    C = np.zeros((num_lvls,num_lvls))
    for p in range(num_partners):
        max_index = len(temps_all[p])
        if T <= temps_all[p][0]: # T is lower than lowest value in table
            for i in range(num_lvls):
                for j in range(num_lvls):
                    C[i][j] = C[i][j] + comp_fractions[p] * C_all[p][0][i][j]
        elif T >= temps_all[p][max_index-1]: # T is higher than highest value in table
            for i in range(num_lvls):
                for j in range(num_lvls):
                    C[i][j] = C[i][j] + comp_fractions[p] * C_all[p][max_index-1][i][j]
        else: # determine temperature entries needed to interpolate
            t = 1 # T index of upper limit
            while temps_all[p][t] < T:
                t = t+1
            t_frac = (temps_all[p][t] - T)/(temps_all[p][t] - temps_all[p][t-1])
            for i in range(num_lvls):
                for j in range(num_lvls):
                    interpol = (1-t_frac) * C_all[p][t][i][j] + t_frac * C_all[p][t-1][i][j]
                    C[i][j] = C[i][j] + comp_fractions[p] * interpol

    return C

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

# def calc_integral1(s, pdf, X_H2_ss, ds):
#     integ1 = 0.0
#     for i in range(0, 100):
#         integ1 += np.exp(s[i]) * pdf[i] * X_H2_ss[i] * ds
#     return integ1

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

# def calc_integral2(s, pdf, X_CO, ds):
#     integ2 = 0.0
#     for i in range(0, 100):
#         integ2 += np.exp(s[i]) * pdf[i] * X_CO[i] * ds
#     return integ2

def calc_line_profile(nu, c_s_CO, c_cgs):
    delta_nu = c_s_CO*nu/c_cgs #"Hz"
    return 1/delta_nu #Hz-1

def tau_LVG(N, nu, lambda_jeans, n_i, n_j, B_ij, B_ji, c_s_CO, c_cgs, h_ev):
    # units: eV*s * Hz * cm * 1/cm3 * cm2/(eV*s) * 1/(Hz) = none
    delta_nu = c_s_CO*nu/c_cgs #"Hz"
    return h_ev*nu*lambda_jeans*N*((n_j*B_ji)-(n_i*B_ij)) / (4*np.pi*delta_nu)
    
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

def calc_integrated_emissivity(N, nu, n_i, A_ij, h_ev):
    #units: eV*s * Hz * cm-3 * s-1 * Hz-1 = eV/cm3
    j_nu = h_ev * nu * N * n_i * A_ij / (4*np.pi)
    return j_nu

# def calc_j_nu_bar(beta_nu, j_nu, pdf, ds):
#     integral3 = 0.0
#     for i in range(0, 100):
#         integral3 += beta_nu[i] * j_nu[i] * pdf[i] * ds 
#     return integral3

def calc_temp_for_rad_field(nu, I_nu, h_ev, c_cgs, K_b_ev):
    foo = ((2*h_ev*(nu**3))/((c_cgs**2)*I_nu)) + 1
    #units: eV*s * Hz3 * cm-2*s2 * eV-1*cm2 = Hz3 * s3
    temp_for_rad_field = h_ev * nu / (np.log(foo) * K_b_ev) #units: eV*s * Hz * eV-1*K = Hz * s * K
    return temp_for_rad_field

def inside_loop(mach_no, n_H_mean, metal, G_o, T_mean, T, G, m_p, K_b, K_b_ev, c_cgs, h_ev, T_bg, eV,
                s, n_H, pdf, lambda_jeans, X_CO, n_CO, sigma_s, s_bar, smin, smax, ds,
                c_s_CO, freq, E, A, B, C_all, C, temps_all, ni, Z, tau_nu, beta_nu, B_nu, j_nu, phi_ij, nu, n_i, n_j ):
    
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
    tau_nu = tau_LVG(N, nu, lambda_jeans, n_i, n_j, B_ij, B_ji, c_s_CO, c_cgs, h_ev)
    beta_nu = beta_LVG(tau_nu)
    j_nu = calc_integrated_emissivity(N, nu, n_i, A_ij, h_ev)
    j_nu_bar = np.sum(beta_nu*j_nu*pdf*ds)
    # temp_nu = pynbody.array.SimArray(calc_temp_for_rad_field(nu, I_nu), "K")
        
    return X_H2_bar, X_CO_bar, j_nu_bar


if __name__=='__main__':
    path = "bulk1/data_2/hydro_59/output/"
    data = pynbody.load(path + "output_00050")
    aexp = data.properties['a']
    data.physical_units()
    r_vir = center_and_r_vir(data, aexp, path)
    remove_bulk_velocity(data)
    r_e = 0.1 * r_vir
    sph_5 = pynbody.filt.Sphere(radius = '5.0 kpc') # %(r_e*1.4))
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
    c_cgs = pynbody.array.SimArray(2.99792458e10, "cm s**-1")
    h_ev = pynbody.array.SimArray(4.13566770e-15, "eV s")
    mH_cgs = pynbody.array.SimArray(1.6737236e-24, 'g')  # hydrogen mass
    T_bg = pynbody.array.SimArray(2.73, "K")
    eV = pynbody.array.SimArray(6.241509e18, "J")

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
    mu, num_lvls, E, g, freq, A, B, C_all, num_partners, temps_all, num_temps, n_coll = read_file('CO.txt',c_cgs, h_ev, K_b_ev)
    comp_fractions, abundance = load_species_info('CO')
    # calc netto collision coeffs
    C = calc_total_C(T_mean, comp_fractions)
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
    j_nu = np.zeros(100)
    u = 1 #upper level
    l = 0 #lower level
    nu = freq[u][l]
    n_i = ni[u]
    n_j = ni[l]
    A_ij = A[u][l]
    B_ij = B[u][l]
    B_ji = B[l][u]
    phi_ij = calc_line_profile(nu, c_s_CO, c_cgs)

    X_H2_bar = np.zeros(len(n_H_mean_arr))
    X_CO_bar = np.zeros(len(n_H_mean_arr))
    j_nu_bar = np.zeros(len(n_H_mean_arr))

    for m in range(0, len(mach_no_arr)):
        mach_no = mach_no_arr[m]
        n_H_mean = n_H_mean_arr[m]
        metal = metal_arr[m]
        T = temperature[m]
        X_H2_bar[m], X_CO_bar[m], j_nu_bar[m] = inside_loop(mach_no, n_H_mean, metal, G_o, T_mean, T, 
                                                            G, m_p, K_b, K_b_ev, c_cgs, h_ev, T_bg, eV,
                                                            s, n_H, pdf, lambda_jeans, X_CO, n_CO, 
                                                            sigma_s, s_bar, smin, smax, ds, 
                                                            c_s_CO, freq, E, A, B, C_all, C, temps_all, ni, Z,
                                                            tau_nu, beta_nu, B_nu, j_nu, phi_ij, nu, n_i, n_j )
    #radiative_transfer/outputs_RT/1.3/
    np.save('radiative_transfer/outputs_RT/1.3/X_H2_bar_1.3.npy', X_H2_bar)
    np.save('radiative_transfer/outputs_RT/1.3/X_CO_bar_1.3.npy', X_CO_bar)
    np.save('radiative_transfer/outputs_RT/1.3/mach_no_arr_1.3.npy', mach_no_arr)
    np.save('radiative_transfer/outputs_RT/1.3/n_H_mean_arr_1.3.npy', n_H_mean_arr)
    np.save('radiative_transfer/outputs_RT/1.3/metal_arr_1.3.npy', metal_arr)
    np.save('radiative_transfer/outputs_RT/1.3/T_1.3.npy', temperature)
    np.save('radiative_transfer/outputs_RT/1.3/j_nu_bar_1.3.npy', j_nu_bar)