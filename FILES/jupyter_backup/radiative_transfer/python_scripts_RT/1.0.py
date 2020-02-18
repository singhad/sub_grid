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
        E.append(float(words[1]) * 100.*c_cgs*h_ev) # cm^-1 -> eV
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

''' calculate the partition function
    PARAMS:
      T = temperature (K)
      num_lvls = number of energy levels
      g = statistical weight of each level
      E = energy of each level (eV)'''
def partion_function(T, num_lvls, g, E):
    Z=0.0
    for i in range(0,num_lvls):
        Z = Z + g[i]*np.exp(-E[i]/(K_b_ev*T))
    return Z

''' calculate LTE occupation numbers with partition function method
    PARAMS:
      T = temperature (K)
      num_lvls = number of energy levels
      g = statistical weight of each level
      E = energy of each level (eV)
    RETURN:
      ni = level populations '''
def calc_lvlpops_partion(T, num_lvls, g, E):
    ni = []
    Z = partion_function(T, num_lvls, g, E)
    for i in range(0, num_lvls):
        ni.append(g[i]*np.exp(-E[i]/(K_b_ev*T)) / Z)
    return ni, Z

''' Optical depth in LVG approximation for line ij '''
def tau_LVG(N, nu_ij, lambda_jeans, n_i, n_j, B_ij, B_ji, c_s):
    # units: eV*s * 1/s * cm * 1/cm3 * cm2/(eV*s) * 1/(1/s) = none
    grad_nu = c_s*nu_ij/c_cgs #units: 1/s
    return h_ev*nu_ij*lambda_jeans*((N*n_j*B_ji)-(N*n_j*B_ij)) / (4*np.pi*grad_nu)

''' Escape probability in LVG approximation for line ij
    PARAMS:
      tau = optical depth'''
def beta_LVG(tau):
    if tau < 0.01:
        return 1. - tau/2.
    elif tau > 100.:
        return 1./tau
    else:
        return (1.0 - np.exp(-tau)) / tau

            ''' total emissivity of transition i-j in W/m3/sr
                J = integral j dv = integral cst phi(v) dv = cst integral phi dv = cst * 1
                Remark that v_ij is constant! '''
            def integrated_emissivity(nu_ij, x_i, n_molec, A_ij):
                # units: J*s * Hz * sr-1 * m-3 * s-1 = J/s/m3/sr
                return h_si*nu_ij/(4*np.pi) * x_i * n_molec * A_ij


            ''' total extinction of transition i-j in 1/m/sr/s
                B in m2/(eV*s) '''
            def integrated_extinction(nu_ij, x_i, x_j, n_molec, B_ij, B_ji):
                # units: eV*s * Hz * sr-1 * m2/(eV*s) * m-3 = 1/m/sr/s
                return h_ev*nu_ij/(4*np.pi) * (x_j*B_ji - x_i*B_ij) * n_molec


            ''' Source function of transision i-j for optically thick media S=j/a in eV/s/Hz/m2 '''
            def source_function_thick(nu_ij, x_i, x_j, n_molec, A_ij, B_ij, B_ji):
                j = integrated_emissivity(nu_ij, x_i, n_molec, A_ij) / eV
                a = integrated_extinction(nu_ij, x_i, x_j, n_molec, B_ij, B_ji)
                # units: (eV/s/m3/sr) / (1/m/sr/s) = eV/m2 = eV/s/Hz/m2
                return j/a

def B_nu_ev(nu, T):
    if nu==0.:
        return 0.
    else:
        x = h_ev*nu/(K_b_ev*T) #units: none
        #units: eV*s * Hz3 / (cm2/s2) = eV * s3 * Hz3 * cm-2 = eV/s/Hz/cm2 = eV/cm2
        return 2.0*h_ev*(nu**3)/(c**2) / (np.exp(x)-1.0)

def LTE_solver(I_nu_bg, B_nu, tau):
    return ((I_nu_bg*np.exp(-tau)) + (B_nu*tau))

for i in range(0, num_lvls):
    for j in range(0, num_lvls):
            if freq[i][j] != 0.0:
                tau_RT[i][j] = tau_LVG(N, freq[i][j], lambda_jeans, ni[i], ni[j], B[i][j], B[j][i], c_s)
                beta_RT[i][j] = beta_LVG(tau_ji[i][j])
                B_nu[i][j] = B_nu_ev(freq[i][j], T_mean)
                I_nu_bg[i][j] = B_nu_ev(freq[i][j], T_bg)
                I_nu[i][j] = LTE_solver(I_nu_bg[i][j], B_nu[i][j], tau_ji[i][j])

if __name__=='__main__':
    n_H_arr = pynbody.array.SimArray(np.load('outputs/sub_grid/n_H.npy'), "cm**-3")
    n_H2_arr = pynbody.array.SimArray(np.load('outputs/sub_grid/n_H2.npy'), "cm**-3")
    X_H2_arr = pynbody.array.SimArray(np.load('outputs/sub_grid/X_H2.npy'), "1")
    n_LW_arr = np.load('outputs/sub_grid/n_LW.npy')
    n_LW_ss_arr = np.load('outputs/sub_grid/n_LW_ss.npy')
    X_H2_ss_arr = pynbody.array.SimArray(np.load('outputs/sub_grid/X_H2_ss.npy'), "1")
    n_H2_ss_arr = pynbody.array.SimArray(np.load('outputs/sub_grid/n_H2_ss.npy'), "cm**-3")
    X_CO_arr = pynbody.array.SimArray(np.load('outputs/sub_grid/X_CO.npy'), "1")
    n_CO_arr = pynbody.array.SimArray(np.load('outputs/sub_grid/n_CO.npy'), "cm**-3")
    pdf_arr = np.load('outputs/sub_grid/pdf.npy')
    lambda_jeans_arr = pynbody.array.SimArray(np.load('outputs/sub_grid/lambda_jeans.npy'), "cm")

    #Defining all the constants used in this whole program

    m_p = pynbody.array.SimArray(1.672621777e-24, "g")
    K_b = pynbody.array.SimArray(1.38064852e-16, "cm**2 g s**-2 K**-1")
    K_b_ev = pynbody.array.SimArray(8.617e-5, "eV K**-1")
    G = pynbody.array.SimArray(6.67259e-8, "cm**3 g**-1 s**-2")
    T_mean = pynbody.array.SimArray(10., "K")
    mach_no = pynbody.array.SimArray(10., "1")
    metallicity = 0.02/0.02
    G_o = 1
    n_H_mean = pynbody.array.SimArray(100., "cm**-3")
    c_si = pynbody.array.SimArray(299792458, "m s**-1")
    c_cgs = pynbody.array.SimArray(29979245800, "cm s**-1")
    h_ev = pynbody.array.SimArray(4.13566770e-15, "eV s")
    h_si = pynbody.array.SimArray(6.626e-34, "J s")
    c_s = np.sqrt(K_b * T_mean/m_p)
    T_bg = pynbody.array.SimArray(2.73, "K")
    eV = pynbody.array.SimArray(6.241509e18, "J")
    #pynbody.array.SimArray(, "")


    mu, num_lvls, E, g, freq, A, B, C_all, num_partners, temps_all, num_temps, n_coll = read_file('CO.txt')

    comp_fractions, abundance = load_species_info('CO')

    # calc netto collision coeffs
    C = calc_total_C(T_mean, comp_fractions)


    ni, Z = calc_lvlpops_partion(T_mean, num_lvls, g, E)
    ni = pynbody.array.SimArray(ni, "1")
    N = n_CO_arr[79]
    Ni = pynbody.array.SimArray(ni*N, "cm**-3")
    n_CO_arr[79]
    lambda_jeans = lambda_jeans_arr[79]
    lambda_jeans

    tau_RT = np.zeros((num_lvls, num_lvls))
    beta_RT = np.zeros((num_lvls, num_lvls))
    B_nu = np.zeros((num_lvls, num_lvls))
    I_nu = np.zeros((num_lvls, num_lvls))
    I_nu_bg = np.zeros((num_lvls, num_lvls))
