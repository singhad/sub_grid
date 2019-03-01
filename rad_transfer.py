import numpy as np
import matplotlib.pyplot as plt
import os
import copy

from constants_lineRT import *

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

def read_file(species):
    filename = get_filename(species)
    f = open(filename, 'r')
    f.readline()
    species = f.readline()
    f.readline()
    mu = float(f.readline()) # molecular weight
    f.readline()
    num_lvls = int(f.readline()) # number of energy levels
    # read energy levels: energy E, statistical weight g
    f.readline()
    E = []
    g = []
    for l in range(num_lvls):
        words = f.readline().split()
        E.append(float(words[1]) * 100.*c_si*h_ev)  # cm^-1 -> eV
        g.append(float(words[2]))
    f.readline()
    num_trans = int(f.readline()) # number of radiative transistions
    # read transistions: upper lvl, lower lvl, A-coefficient, frequency
    f.readline()
    A = np.zeros((num_lvls, num_lvls))
    freq = np.zeros((num_lvls, num_lvls))
    for t in range(num_trans):
        words = f.readline().split()
        i = int(words[1]) - 1
        j = int(words[2]) - 1
        A[i][j] = float(words[3])  # s^-1
        freq[i][j] = float(words[4]) #* 1e9  # GHz -> Hz
        freq[j][i] = freq[i][j]
    # compute B-coefficient via Einstein relations
    # Bij = coeff for stimulated emission, Bji = coeff for extinction (j<i)
    B = np.zeros((num_lvls, num_lvls))
    for i in range(0, num_lvls):
        for j in range(0, i):
            if A[i][j] != 0:
                B[i][j] = A[i][j] * (c_si**2) / \
                    (2*h_ev * (freq[i][j])**3)  # m2/(eV*s)
                B[j][i] = B[i][j] * g[i]/g[j]

    return mu, num_lvls, np.array(E), np.array(g), freq, A, B

def partion_function(T, num_lvls, g, E):
    Z = 0.0
    for i in range(0, num_lvls):
        Z = Z + g[i]*np.exp(-E[i]/(kb_ev*T))
    return Z

def calc_lvl_pops_partion(T, num_lvls, g, E):
    ni = []
    Z = partion_function(T, num_lvls, g, E)
    for i in range(0, num_lvls):
        ni.append(g[i]*np.exp(-E[i]/(kb_ev*T)) / Z)
    return np.array(ni)

def calc_source_func(T, nu_ij):
    #since this is LTE, source function = BBR
    if nu_ij==0.:
        return 0.
    else:
        x = h_si*nu_ij/(kb_si*T)
        # units eV*s * Hz3 / (m2/s2) = eV * s3 * Hz3 * m-2 = eV/s/Hz/m2
        S_nu = 2.0*h_ev*(nu_ij**3)/(c**2) / (np.exp(x)-1.0)
        return S_nu

def calc_extinction(T, n_j, nu_ij, B_ji):
    x = h_si*nu_ij/(kb_si*T)
    return (h_ev/(4*np.pi)) * nu_ij * n_j * B_ji * (1.0 - np.exp(-1*x))

def calc_emissivity(extinction, S):
    return extinction*S

def make_pdf(s, s_bar, sigma_s):
    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_bar)/sigma_s)**2)))
    return pdf

def calc_lambda_jeans(n_H):
    m_p = 1.672621777e-24   # g
    T_mean = 10.            #K
    K_b = 1.38064852e-16    # ergs K-1
    G = 6.67408e-8          # dyne cm^2 g^-2
    lambda_jeans = ((np.sqrt(K_b * T_mean / m_p)) / np.sqrt(4* np.pi * G * n_H * m_p))
    return lambda_jeans

def calc_n_LW(n_H, G_o, lambda_jeans):
    m_p = 1.672621777e-24   # g
    kappa = 1000 * m_p
    rad_field_outside = G_o #in solar units
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    n_LW = rad_field_outside * exp_tau
    return n_LW

def calc_X_H2(n_H, Z, n_LW):
    DC = 1.7e-11
    CC = 2.5e-17            #cm3 s-1
    numerator = DC * n_LW
    denominator = CC * Z * n_H
    X_H2 = 1 / (2 + (numerator/denominator) )
    return X_H2

def calc_stats():
    s = np.zeros(1000)
    pdf = np.zeros(1000)
    lambda_jeans = np.zeros(1000)
    X_H2 = np.zeros(1000)
    n_H = np.zeros(1000)
    n_LW = np.zeros(1000)
    n_H2 = np.zeros(1000)
    mach_no = 5
    sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -7*sigma_s + s_bar
    smax = 7*sigma_s + s_bar
    ds = (smax - smin)/1000
    n_H_mean = 100
    Z = 1
    G_o = 1
    for i in range(0, 1000):
        s[i] = smin + i*ds
    pdf = make_pdf(s, s_bar, sigma_s)
    n_H = n_H_mean * np.exp(s)
    lambda_jeans = calc_lambda_jeans(n_H)
    n_LW = calc_n_LW(n_H, G_o, lambda_jeans)
    X_H2 = calc_X_H2(n_H, Z, n_LW)
    n_H2 = n_H * X_H2

    return n_H, n_H2, X_H2, lambda_jeans

def plotting(L, n_H, ni, source_func, extinction, emissivity, n_H2, X_H2):
    i = 2
    j = 1
    #for m in range(0,1000):
    plt.plot(n_H, L[i][j])
    plt.xlabel('n_H')
    plt.ylabel('L/$L_{sun}$')
    plt.grid(b=True, which='both', axis='both')
    plt.title('n_H vs L/$L_{sun}$')
    plt.savefig(os.path.join('n_HvsL.png'.format()))
    plt.clf()

    plt.plot(np.log10(n_H), L[i][j])
    plt.xlabel('log(n_H)')
    plt.ylabel('L/$L_{sun}$')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs L/$L_{sun}$')
    plt.savefig(os.path.join('log(n_H)vsL.png'.format()))
    plt.clf()

    plt.plot(np.log10(n_H), np.log10(L[i][j]))
    plt.xlabel('log(n_H)')
    plt.ylabel('log(L/$L_{sun}$')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs log(L/$L_{sun}$')
    plt.savefig(os.path.join('log(n_H)vslog(L).png'.format()))
    plt.clf()

    for z in range(0,1000):
        plt.plot(ni[i], L[i][j][z])
    plt.xlabel('ni')
    plt.ylabel('L/$L_{sun}$')
    plt.grid(b=True, which='both', axis='both')
    plt.title('ni vs L/$L_{sun}$')
    plt.savefig(os.path.join('nivsL.png'.format()))
    plt.clf()


if __name__ == '__main__':
    path = 'for RT_LTE'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    T = 10 #temperature (K)
    mu, num_lvls, E, g, freq, A, B = read_file('CO.txt')
    source_func = np.zeros((num_lvls, num_lvls))
    extinction = np.zeros((num_lvls, num_lvls))
    emissivity = np.zeros((num_lvls, num_lvls))
    L = np.zeros((num_lvls, num_lvls, 1000))
    L_sun = 3.828e26
    #calculate LTE occupation numbers with partition function method
    ni = calc_lvl_pops_partion(T, num_lvls, g, E)
    for i in range(0, num_lvls):
        for j in range(0, i):
            source_func[i][j] = calc_source_func(T, freq[i][j])
            extinction[i][j] = calc_extinction(T, ni[j], freq[i][j], B[j][i])
            emissivity[i][j] = calc_emissivity(extinction[i][j], source_func[i][j])
    n_H, n_H2, X_H2, lambda_jeans = calc_stats()
    for i in range(0, num_lvls):
        for j in range(0, i):
            for m in range(0,1000):
                L[i][j][m] = (4*np.pi)* ((lambda_jeans[m])**2) * source_func[i][j] / L_sun

    plotting(L, n_H, ni, source_func, extinction, emissivity, n_H2, X_H2)
