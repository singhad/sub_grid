import numpy as np
import matplotlib.pyplot as plt
import os
import copy
from matplotlib.ticker import FormatStrFormatter as ticks

from constants_lineRT import *
from X_H2_self_shielding import *

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
    mu = float(f.readline())  # molecular weight
    f.readline()
    num_lvls = int(f.readline())  # number of energy levels
    # read energy levels: energy E, statistical weight g
    f.readline()
    E = []
    g = []
    for l in range(num_lvls):
        words = f.readline().split()
        E.append(float(words[1]) * 100.*c_si*h_ev)  # cm^-1 -> eV
        g.append(float(words[2]))
    f.readline()
    num_trans = int(f.readline())  # number of radiative transistions
    # read transistions: upper lvl, lower lvl, A-coefficient, frequency
    f.readline()
    A = np.zeros((num_lvls, num_lvls))
    freq = np.zeros((num_lvls, num_lvls))
    for t in range(num_trans):
        words = f.readline().split()
        i = int(words[1]) - 1
        j = int(words[2]) - 1
        A[i][j] = float(words[3])  # s^-1
        freq[i][j] = float(words[4]) * 1e9  # GHz -> Hz
        freq[j][i] = freq[i][j]
    # compute B-coefficient via Einstein relations
    # Bij = coeff for stimulated emission, Bji = coeff for extinction (j<i)
    B = np.zeros((num_lvls, num_lvls))
    for i in range(0, num_lvls):
        for j in range(0, i):
            if A[i][j] != 0:
                B[i][j] = A[i][j] * (c_cgs**2) / \
                    (2*h_ev * (freq[i][j])**3)  # m2/(eV*s)
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
        num_collis = int(f.readline())
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
        for col in range(num_collis):
            words = f.readline().split()
            i = int(words[1]) - 1
            j = int(words[2]) - 1
            for t in range(num_temps):
                C[t][i][j] = float(words[3+t])  # * 1.e-6 # cm3/s -> m3/s
        # calculate the inverse coefficient via LTE relation
        for i in range(num_lvls):
            for j in range(i):
                for t in range(num_temps):
                    if C[t][i][j] != 0:
                        C[t][j][i] = C[t][i][j] * \
                            np.exp(-(E[i]-E[j])/(kb_ev*temps[t]))*g[i]/g[j]
        # add collision partner data to global array
        C_all.append(C)

    f.close()
    C_all = np.array(C_all)
    temps_all = np.array(temps_all)
    return mu, num_lvls, np.array(E), np.array(g), freq, A, B, C_all, num_partners, temps_all


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


def calc_extinction(T, n_i, n_j, nu_ij, B_ij, B_ji):
    x = h_ev*nu_ij/(kb_ev*T)
    #(h_ev/(4*np.pi)) * nu_ij * n_j * B_ji * (1.0 - np.exp(-1*x))
    extinction = h_ev * nu_ij * (n_j*B_ji - n_i*B_ij) / (4*np.pi)
    return extinction


def calc_emissivity(nu_ij, n_i, A_ij):
    emissivity = h_ev * nu_ij * n_i * A_ij / (4*np.pi)
    return emissivity


def calc_source_func(emissivity, extinction):
    '''#since this is LTE, source function = BBR
    if nu_ij==0.:
        return 0.
    else:
        x = h_ev*nu_ij/(kb_ev*T)
        # units eV*s * Hz3 / (cm2/s2) = eV * s3 * Hz3 * cm-2 = eV/s/Hz/cm2
        S_nu = 2.0*h_ev*(nu_ij**3)/(c_cgs**2) / (np.exp(x)-1.0)
        return S_nu'''
    S_nu = emissivity / extinction
    return S_nu


def make_pdf(s, s_bar, sigma_s):
    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * \
        (np.exp(-0.5*(((s - s_bar)/sigma_s)**2)))
    return pdf


def calc_lambda_jeans(n_H):
    T_mean = 10.  # K
    lambda_jeans = ((np.sqrt(kb_ev * T_mean / m_p)) /
                    np.sqrt(4 * np.pi * G_cgs * n_H * m_p))
    return lambda_jeans


def calc_n_LW(n_H, G_o, lambda_jeans):
    kappa = 1000 * m_p
    rad_field_outside = G_o  # in solar units
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    n_LW = rad_field_outside * exp_tau
    return n_LW


def calc_n_LW_ss(n_H, n_H2, G_o, lambda_jeans):
    m_p = 1.672621777e-24   # g
    kappa = 1000 * m_p
    rad_field_outside = G_o  # in solar units
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    N_H2 = n_H2*lambda_jeans
    term1 = (0.965/((1+(N_H2/5e14))**2))
    term2 = ((0.035/np.sqrt(1+(N_H2/5e14))) *
             np.exp(-1*np.sqrt(1+(N_H2/5e14))/1180))
    S_H2 = term1 + term2
    n_LW_ss = rad_field_outside * exp_tau * S_H2
    return n_LW_ss, S_H2, N_H2


def calc_X_H2(n_H, Z, n_LW):
    DC = 1.7e-11
    CC = 2.5e-17  # cm3 s-1
    numerator = DC * n_LW
    denominator = CC * Z * n_H
    X_H2 = 1 / (2 + (numerator/denominator))
    return X_H2


def calc_stats():
    s = np.zeros(1000)
    pdf = np.zeros(1000)
    lambda_jeans = np.zeros(1000)
    X_H2 = np.zeros(1000)
    n_H = np.zeros(1000)
    n_LW = np.zeros(1000)
    n_LW_ss = np.zeros(1000)
    n_H2 = np.zeros(1000)
    S_H2 = np.zeros(1000)
    N_H2 = np.zeros(1000)
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
        X_H2_a = calc_X_H2(n_H, Z, n_LW)
        n_H2_a = n_H * X_H2_a
        n_LW_1, S_H2_1, N_H2_1 = calc_n_LW_ss(n_H, n_H2_a, G_o, lambda_jeans)
        X_H2_1 = calc_X_H2(n_H, Z, n_LW_1)
        n_H2_1 = n_H * X_H2_1
        n_LW_2, S_H2_2, N_H2_2 = calc_n_LW_ss(n_H, n_H2_1, G_o, lambda_jeans)
        X_H2_2 = calc_X_H2(n_H, Z, n_LW_2)
        n_H2_2 = n_H * X_H2_2
        n_LW_3, S_H2_3, N_H2_3 = calc_n_LW_ss(n_H, n_H2_2, G_o, lambda_jeans)
        X_H2_3 = calc_X_H2(n_H, Z, n_LW_3)
        n_H2_3 = n_H * X_H2_3
        n_LW_4, S_H2_4, N_H2_4 = calc_n_LW_ss(n_H, n_H2_3, G_o, lambda_jeans)
        X_H2_4 = calc_X_H2(n_H, Z, n_LW_4)
        n_H2_4 = n_H * X_H2_4
        n_LW_5, S_H2_5, N_H2_5 = calc_n_LW_ss(n_H, n_H2_4, G_o, lambda_jeans)
        X_H2_5 = calc_X_H2(n_H, Z, n_LW_5)
        n_H2_5 = n_H * X_H2_5
        n_LW_6, S_H2_6, N_H2_6 = calc_n_LW_ss(n_H, n_H2_5, G_o, lambda_jeans)
        X_H2_6 = calc_X_H2(n_H, Z, n_LW_6)
        n_H2_6 = n_H * X_H2_6
        n_LW_7, S_H2_7, N_H2_7 = calc_n_LW_ss(n_H, n_H2_6, G_o, lambda_jeans)
        X_H2_7 = calc_X_H2(n_H, Z, n_LW_7)
        n_H2_7 = n_H * X_H2_7
        n_LW_8, S_H2_8, N_H2_8 = calc_n_LW_ss(n_H, n_H2_7, G_o, lambda_jeans)
        X_H2_8 = calc_X_H2(n_H, Z, n_LW_8)
        n_H2_8 = n_H * X_H2_8
        n_LW_9, S_H2_9, N_H2_9 = calc_n_LW_ss(n_H, n_H2_8, G_o, lambda_jeans)
        X_H2_9 = calc_X_H2(n_H, Z, n_LW_9)
        n_H2_9 = n_H * X_H2_9
        n_LW_10, S_H2_10, N_H2_10 = calc_n_LW_ss(
            n_H, n_H2_9, G_o, lambda_jeans)
        X_H2_10 = calc_X_H2(n_H, Z, n_LW_10)
        n_H2_10 = n_H * X_H2_10
        n_LW_ss, S_H2, N_H2 = calc_n_LW_ss(n_H, n_H2_10, G_o, lambda_jeans)
        X_H2 = calc_X_H2(n_H, Z, n_LW_ss)
        n_H2 = n_H * X_H2

    return n_H, n_H2, X_H2, lambda_jeans


def plotting(L, n_H, ni, source_func, extinction, emissivity, n_H2, X_H2, freq):
    fig, ax = plt.subplots()
    i = 2
    j = 1
    '''m = 1

    for u in range(0, num_lvls):
        for l in range(0, u):
            if u-l==1:
                plt.scatter(n, emissivity[u][l])
    plt.xlabel('ni')
    plt.ylabel('Emissivity')
    plt.grid(b=True, which='both', axis='both')
    ax.yaxis.set_major_formatter(ticks('%.2f'))
    plt.title('ni vs Emissivity')
    plt.savefig(os.path.join('nivsEmissivity'.format()), dpi= 'figure', bbox_inches= 'tight')
    plt.clf()'''

    # for m in range(0,1000):
    plt.plot(n_H, L[i][j])
    plt.xlabel('$n_H$')
    plt.ylabel('L/$L_{sun}$')
    plt.grid(b=True, which='both', axis='both')
    ax.yaxis.set_major_formatter(ticks('%.2f'))
    plt.title('n_H vs (L/$L_{sun}$)')
    plt.savefig(os.path.join('n_HvsL.png'.format()),
                dpi='figure', bbox_inches='tight')
    plt.clf()

    plt.plot(np.log10(n_H), L[i][j])
    plt.xlabel('log($n_H$)')
    plt.ylabel('(L/$L_{sun}$)')
    plt.grid(b=True, which='both', axis='both')
    ax.yaxis.set_major_formatter(ticks('%.2f'))
    plt.title('log(n_H) vs (L/$L_{sun}$)')
    plt.savefig(os.path.join('log(n_H)vsL.png'.format()),
                dpi='figure', bbox_inches='tight')
    plt.clf()

    for u in range(0, num_lvls):
        for l in range(0, u):
            if u-l == 1:
                plt.scatter(np.log10(freq[u][l]), source_func[u][l])
    plt.xlabel('freq')
    plt.ylabel('$S_{nu}$')
    plt.grid(b=True, which='both', axis='both')
    ax.yaxis.set_major_formatter(ticks('%.2f'))
    plt.title('freq vs $S_{nu}$')
    plt.savefig(os.path.join('freqvsS_nu'.format()),
                dpi='figure', bbox_inches='tight')
    plt.clf()

    plt.plot(np.log10(n_H), np.log10(L[i][j]))
    plt.xlabel('log(n_H)')
    plt.ylabel('log(L/$L_{sun}$)')
    plt.grid(b=True, which='both', axis='both')
    ax.yaxis.set_major_formatter(ticks('%.1f'))
    plt.title('log(n_H) vs log(L/$L_{sun}$)')
    plt.savefig(os.path.join('log(n_H)vslog(L).png'.format()))
    plt.clf()


if __name__ == '__main__':
    path = 'for RT_LTE'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)
    T = 10  # temperature (K)
    mu, num_lvls, E, g, freq, A, B, C, num_partners, temps_all = read_file('CO.txt')
    source_func = np.zeros((num_lvls, num_lvls))
    extinction = np.zeros((num_lvls, num_lvls))
    emissivity = np.zeros((num_lvls, num_lvls))
    L = np.zeros((num_lvls, num_lvls, 1000))
    # calculate LTE occupation numbers with partition function method
    ni = calc_lvl_pops_partion(T, num_lvls, g, E)
    for i in range(0, num_lvls):
        for j in range(0, i):
            extinction[i][j] = calc_extinction(
                T, ni[i], ni[j], freq[i][j], B[i][j], B[j][i])
            emissivity[i][j] = calc_emissivity(freq[i][j], ni[i], A[i][j])
            source_func[i][j] = calc_source_func(
                emissivity[i][j], extinction[i][j])
    n_H, n_H2, X_H2, lambda_jeans = calc_stats()
    for i in range(0, num_lvls):
        for j in range(0, i):
            for m in range(0, 1000):
                L[i][j][m] = (4*np.pi) * ((lambda_jeans[m]) **
                                          2) * source_func[i][j] / L_sun

    plotting(L, n_H, ni, source_func, extinction, emissivity, n_H2, X_H2, freq)
