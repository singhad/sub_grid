
# coding: utf-8

# In[ ]:


# File for testing different RT codes before applying them to python scripts


# In[1]:


import sys
import numpy as np
np.set_printoptions(threshold=sys.maxsize)
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import pynbody
from michaels_functions import center_and_r_vir, remove_bulk_velocity
from matplotlib.colors import LogNorm
from matplotlib.pyplot import figure


# In[2]:


n_H_arr = np.load('outputs/sub_grid/n_H.npy')
n_H2_arr = np.load('outputs/sub_grid/n_H2.npy')
X_H2_arr = np.load('outputs/sub_grid/X_H2.npy')
n_LW_arr = np.load('outputs/sub_grid/n_LW.npy')
n_LW_ss_arr = np.load('outputs/sub_grid/n_LW_ss.npy')
X_H2_ss_arr = np.load('outputs/sub_grid/X_H2_ss.npy')
n_H2_ss_arr = np.load('outputs/sub_grid/n_H2_ss.npy')
X_CO_arr = np.load('outputs/sub_grid/X_CO.npy')
n_CO_arr = np.load('outputs/sub_grid/n_CO.npy')
pdf_arr = np.load('outputs/sub_grid/pdf.npy')
lambda_jeans_arr = np.load('outputs/sub_grid/lambda_jeans.npy')


# In[5]:


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
        E.append(float(words[1]) * 100.*c_si*h_ev) # cm^-1 -> eV
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
                B[i][j] = A[i][j] * (c_si**2) / (2*h_ev * (freq[i][j])**3) # m2/(eV*s)
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
                C[t][i][j] = float(words[3+t]) * 1.e-6 # cm3/s -> m3/s
        
        # calculate the inverse coefficient via LTE relation
        for i in range(num_lvls):
            for j in range(i):
                for t in range(num_temps):
                    if C[t][i][j] != 0:
                        C[t][j][i] = C[t][i][j] * np.exp(-(E[i]-E[j])/(kb_ev*temps[t]))*g[i]/g[j]
        
        # add collision partner data to global array
        C_all.append(C)

    f.close()
    C_all = pynbody.array.SimArray(np.array(C_all), "m**3 s**-1")
    temps_all = pynbody.array.SimArray(np.array(temps_all), "K")
    E = pynbody.array.SimArray(np.array(E), "eV")
    g = pynbody.array.SimArray(np.array(g), "1")
    freq = pynbody.array.SimArray(np.array(freq), "Hz")
    A = pynbody.array.SimArray(np.array(A), "s**-1")
    B = pynbody.array.SimArray(np.array(B), "m**2 eV**-1 s**-1")
    return mu, num_lvls, E, g, freq, A, B, C_all, num_partners, temps_all, num_temps, num_coll

#pynbody.array.SimArray(, "")


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
                    interpol = (1-t_frac) * C_all[p][t][i][j] +                                    t_frac * C_all[p][t-1][i][j]
                    C[i][j] = C[i][j] + comp_fractions[p] * interpol
                    
    return C


# In[7]:


#------background spectra---------------------------------------------------------------------------
''' Planck spectrum for a black body
    PARAMS:
      nu = freq in Hz
      T = temperature in K
    RETURN:
      spectral radiance in W/sr/m2/Hz '''
def B_nu_si(nu, T):
    x = h_si*nu/(kb_si*T)
    # units J/Hz * Hz3 * s2/m2 = J/Hz * Hz * m-2 = J/s /Hz /m2 = W/(m2*Hz)
    return 2.0*h_si*(nu**3)/(c**2) / (np.exp(x)-1.0)

def B_nu_ev(nu, T):
    if nu==0.:
        return 0.
    else:
        x = h_si*nu/(kb_si*T)
        # units eV*s * Hz3 / (m2/s2) = eV * s3 * Hz3 * m-2 = eV/s/Hz/m2
        return 2.0*h_ev*(nu**3)/(c**2) / (np.exp(x)-1.0)
    
    
#-----Definitions------------------------------------------------------------------------------------

''' total emissivity of transition i-j in W/m3/sr
    J = integral j dv = integral cst phi(v) dv = cst integral phi dv = cst * 1
    Remark that v_ij is constant! '''
def integrated_emissivity(nu_ij, n_i, n_molec, A_ij):
    # units: J*s * Hz * sr-1 * m-3 * s-1 = J/s/m3/sr
    return h_si*nu_ij/(4*np.pi) * n_i * n_molec * A_ij


''' total extinction of transition i-j in 1/m/sr/s
    B in m2/(eV*s) '''
def integrated_extinction(nu_ij, n_i, n_j, n_molec, B_ij, B_ji):
    # units: eV*s * Hz * sr-1 * m2/(eV*s) * m-3 = 1/m/sr/s
    return h_ev*nu_ij/(4*np.pi) * (n_j*B_ji - n_i*B_ij) * n_molec


''' Source function of transision i-j for optically thick media S=j/a in eV/s/Hz/m2 '''
def source_function_thick(nu_ij, n_i, n_j, n_molec, A_ij, B_ij, B_ji):
    j = integrated_emissivity(nu_ij, n_i, n_molec, A_ij) / eV  #converting from J to eV
    a = integrated_extinction(nu_ij, n_i, n_j, n_molec, B_ij, B_ji)
    # units: (eV/s/m3/sr) / (1/m/sr/s) = eV/m2 = eV/s/Hz/m2
    return j/a

''' Thermal line profile
    PARAMS:
      freq = frequency (Hz)
      centerFreq = rest frequency of the transition (Hz)
      T = temperature (K)
      moleculeMass = the mass of the particle (atomic units)
      v_los (optional) = the line of sight velocity to calc the doppler shift '''
def thermal_profile(freq, freq_center, T, moleculeMass, vlos=0.):
    if (freq==0.) or (freq_center==0.):
        profile = 0.
    else:
        new_center = freq_center * (1. + v_los/c_si)
        # thermal line width in units: sqrt(J/K *K /kg) = sqrt(kg m2/s2 /kg) =m/s
        a = np.sqrt(2.*kb_si*T/(moleculeMass*MH))
        f = c_si/(a*new_center*np.sqrt(np.pi)) 
        # units: m/s * Hz / (m/s *Hz) * m/s /(m/s * Hz) = 1/Hz
        profile = f * np.exp(-(c_si*abs(freq-new_center)/(a*new_center))**2)
    return profile


''' Integrated line profile over [freq_a, freq_b]
    This is basically the integration of a Gaussian distribution
    C sqrt(pi) int_p^q exp(-(C*x)**2) dx = 0.5 * (erf(q*C) - erf(p*C))
    with C = 1/((sqrt(2)*sigma) and sigma = therm_a * freq_center/(sqrt(2)*c)'''
def integrate_thermal_profile(freq_a, freq_b, freq_center, T, moleculeMass, v_los=0.):
    # doppler shift
    new_center = freq_center * (1. + v_los/c_si)
    # shift integration limits so line is centered on around 0
    a = freq_a - new_center
    b = freq_b - new_center
    # calc the thermal width of the line
    therm_a = np.sqrt(2.*kb_si*T/(moleculeMass*MH))
    # calc the integral
    C = c_si / (therm_a * freq_center)
    return 0.5 * (erf(b*C) - erf(a*C))

def test_integration():
    freqs = np.linspace(1.99999e9,2.00001e9,10)
    freq_center = 2.e9
    T = 10.
    mu = 29.
    therm_a = np.sqrt(2.*kb_si*T/(mu*MH))
    tot = 0.
    for f in range(len(freqs)-1):
        part = integrate_thermal_profile(freqs[f], freqs[f+1], freq_center, T, mu, 0.)
        print part
        tot = tot + part
    print tot
    
#---------------------------------------------------------------------------------------------------
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
    return ni


# In[8]:


''' Full non-LTE, with radiation field:
    Calculate non-LTE occupation numbers by solving the system of coupled balance equations
    PARAMS:
      num_lvls = number of energy levels
      A = Einstein A coeff matrix (1/s)
      B = Einstein B coeff matrix (m2/(eV*s))
      C = netto collision coeff matrix (m3/s)
      n_coll = total number density of collision partners (1/m3)
      rad_field = incomming radiation field for each transition I[i][j] (eV/s/m2/Hz = eV/m2)
    RETURN:
      sol = level populations '''
def calc_lvlpops_nonLTE(num_lvls, A, B, C, n_coll, rad_field):
    # solve  M*n = 0
    # n = [n1,n2,...ni,...,nn]

    # fill matrix M
    M = np.zeros((num_lvls,num_lvls))
    for a in range(0, num_lvls):
        for b in range(0,num_lvls):
            M_ab = 0
            # upper triangle
            if b>a:
                M_ab = A[b][a] + B[b][a]*rad_field[b][a] + C[b][a]*n_coll #1/s
            # diagonal
            elif a==b:
                for j in range(0, a):
                    M_ab = M_ab - A[a][j] - (B[a][j]*rad_field[a][j]) - C[a][j]*n_coll
                for j in range(a+1, num_lvls):
                    M_ab = M_ab - B[a][j]*rad_field[j][a] - C[a][j]*n_coll
            # lower triangle
            else:
                M_ab = B[b][a]*rad_field[a][b] + C[b][a]*n_coll
            M[a][b] = M_ab

    # solve M*n=0 with svd: M = U*S*V.T
    U, S, Vt = np.linalg.svd(M)
    # In S the smallest singular value is given last. Check if it sufficiently small
    if S[num_lvls-1] > 1.0e-4*S[num_lvls-2]:
        print 'WARNING: unreliable solution:', S 
    sol = Vt[num_lvls-1]

    # sol can be multiplied by a constant factor: normalise to sum(sol_i) = 1
    # (assumes all values are either positive or neg)
    norm=np.sum(sol)
    sol=sol/norm
    return sol

''' Self-consistently solve for the level populations with the general method
    in the optical thick limit!
    PARAMS:
      n_coll = total number density of collision partners (1/m3)
      comp_fractions = fractions of n_coll for each collision partner (from species file)
      n_molec = number density of the molecule (1/m3)
      T = temperature (K)
      rad_field_bg =  background radiation field (CMB)
      num_iter = number of iterations for determining the lvl pops'''
def solve_lvlpops_nonLTE(n_coll, comp_fractions, n_molec, T, rad_field_bg, num_iter=10):

    # calc netto collision coeffs
    C = calc_total_C(T, comp_fractions)

    # set initial radiation field
    source = np.zeros((num_lvls,num_lvls))
    #J = emis * 4.*np.pi #* dx # eV/m2
    J = rad_field_bg

    for it in range(num_iter):
        lvl_pops = calc_lvlpops_nonLTE(num_lvls, A, B, C, n_coll, J)
        # calculate the emissivity
        for i in range(num_lvls):
            for j in range(i):
                if freq[i][j] != 0.:
                    source[i][j] = source_function_thick(freq[i][j],
                                             lvl_pops[i], lvl_pops[j],
                                             n_molec, A[i][j],
                                             B[i][j], B[j][i])
        J = source
        #J = J + rad_field_bg

    return lvl_pops



# In[9]:


#-------LVG solver----------------------------------------------------------------------------------
""" LVG + EscProb method 
    Calculate non-LTE occupation numbers by solving the system of coupled balance equations with
    the LVG approximation
    PARAMS:
      num_lvls = number of energy levels
      A = Einstein A coeff matrix (1/s)
      B = Einstein B coeff matrix (m2/(eV*s))
      C = netto collision coeff matrix (m3/s)
      n_coll = total number density of collision partners (1/m3)
      rad_field_bg = background radiation field for each transition I[i][j] (eV/s/m2/Hz = eV/m2)
      beta =  escape probability for each transition (/)
    RETURN:
      sol = level populations"""
def calc_lvlpops_LVG(num_lvls, A, B, C, n_coll, rad_field_bg, beta):
    # solve  M*n = 0
    # n = [n1,n2,...ni,...,nn]

    # fill matrix M
    M = np.zeros((num_lvls,num_lvls))
    for a in range(0, num_lvls):
        for b in range(0,num_lvls):
            M_ab = 0
            # upper triangle
            if b>a:
                M_ab = A[b][a]*beta[b][a] + B[b][a]*rad_field_bg[b][a]*beta[b][a] + C[b][a]*n_coll
            # diagonal
            elif a==b:
                for j in range(0, a):
                    M_ab = M_ab - A[a][j]*beta[a][j] - B[a][j]*rad_field_bg[a][j]*beta[a][j]                           - C[a][j]*n_coll
                for j in range(a+1, num_lvls):
                    M_ab = M_ab - B[a][j]*rad_field_bg[j][a]*beta[j][a] - C[a][j]*n_coll
            # lower triangle
            else:
                M_ab = B[b][a]*rad_field_bg[a][b]*beta[a][b] + C[b][a]*n_coll
            M[a][b] = M_ab

    # solve M*n=0 with svd: M = U*S*V.T
    U, S, Vt = np.linalg.svd(M)
    # In S the smallest singular value is given last. Check if it sufficiently small
    if S[num_lvls-1] > 1.0e-4*S[num_lvls-2]:
        print 'WARNING: unreliable solution', S
    sol = Vt[num_lvls-1]

    # sol can be multiplied by a constant factor: normalise to sum(sol_i) = 1
    # (assumes all values are either positive or neg)
    norm=np.sum(sol)
    sol=sol/norm
    return sol

''' Self-consistently solve for the level populations with LVG approximation 
    PARAMS:
      line_data = lineData object with info about the transitions
      n_coll = total number density of collision partners (1/m3)
      comp_fracs = fractions of n_coll for each collision partner
      n_molec = number density of the molecule (1/m3)
      T = temperature (K)
      grad_v = velocity gradient (s-1)
      rad_field_bg =  background radiation field (CMB)
      num_iter = maximum number of iterations for determining the lvl pops'''
def solve_lvlpops_LVG(n_coll, comp_fracs, n_molec, T, grad_v, rad_field_bg, num_iter):
    # calc netto collision coeffs
    C = calc_total_C(T, comp_fracs)

    # initialize x_i to LTE
    lvl_pops = calc_lvlpops_partion(T, num_lvls, g, E)

    # calculate tau and beta for LVG
    tau, beta = LVG_coeffs_all_lines(n_molec, lvl_pops, B, grad_v)
    T_ex_prev = 0.0
    T_ex_curr = calc_T_ex(1, 0, lvl_pops, g, E)
    it=0
    convergence = abs(T_ex_curr-T_ex_prev)/(T_ex_curr+T_ex_prev)/2.

    # update level pops with LVG and iterate
    while (convergence>0.01) and (it<num_iter):
        it=it+1
        lvl_pops = calc_lvlpops_LVG(num_lvls, A, B, C, n_coll, rad_field_bg, beta)
        tau, beta = LVG_coeffs_all_lines(n_molec, lvl_pops, B, grad_v)
        T_ex_prev = T_ex_curr
        T_ex_curr = calc_T_ex(1, 0, lvl_pops, g, E)
        convergence = abs(T_ex_curr-T_ex_prev)/(T_ex_curr+T_ex_prev)/2.
        #print 'iteration {} -> convergence {:1.3}%'.format(it, convergence)

    return lvl_pops, beta

''' Calculate tau and beta in LVG approx for all lines
    PARAMS:
      n_molec = number density of the molecule (1/m3)
      lvl_pops = level populations
      B = Einstein B coeff matrix (m2/(eV*s))
      grad_v = velocity gradient (1/s) '''
def LVG_coeffs_all_lines(n_molec, lvl_pops, B, grad_v):
    ni, nj = B.shape
    tau = np.ones((ni,nj))
    beta = np.ones((ni,nj))
    for i in range(ni):
        for j in range(i):
            if B[i][j]==0:
                tau[i][j] = 1.0 # arbitrary
                beta[i][j] = 1.0 #0.0
            else:
                tau[i][j] = tau_LVG(n_molec, grad_v, lvl_pops[i], lvl_pops[j], B[i][j], B[j][i])
                beta[i][j] = beta_LVG(tau[i][j])
                #print 'i {} j {} tau={} beta={}'.format(i,j,tau[i][j],beta[i][j]) 
    return tau, beta

''' Optical depth in LVG approximation for line ij '''
def tau_LVG(n_molec, grad_v, x_i, x_j, B_ij, B_ji):
    # units: m/s * eV*s * 1/m3 * m2/(eV*s) / (1/s) = none
    return c*h_ev/(4.*np.pi) * n_molec * (x_j*B_ji - x_i*B_ij)  / (1.064*grad_v)

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

''' Determine the density for which tau LVG is 1 using B coeffs '''
def critical_density_B(grad_v, x_i, x_j, B_ij, B_ji, abundance):
    n_molec_crit = (4.*np.pi*1.064*grad_v)/((x_j*B_ji - x_i*B_ij)*c*h)
    print 'n_crit', n_molec_crit/ 1e6 / abundance, 'H2/cc'
    return n_molec_crit / abundance * (2*MH) # kg/m3

''' Determine the density for which tau LVG is 1 using A coeff '''
def critical_density_A(grad_v, x_i, x_j, A_ij, freq, g_i, g_j, abundance):
    n_molec_crit = (8.*np.pi*(freq**3)*1.064*grad_v)/((x_j*(g_i/g_j) - x_i)*A_ij*(c**3))
    print 'n_crit', n_molec_crit/ 1e6 / abundance, 'H2/cc'
    print 'args rho_crit:', grad_v, x_i, x_j, A_ij, freq, g_i, g_j, abundance
    return n_molec_crit / abundance * (2*MH) # kg/m3

''' TODO Estimate the critical density in function of grad_v in the simulation '''

#---------------------------------------------------------------------------------------------------
''' Calculate the excitation temperature. This is the LTE temperature corresponding to the given x_i
    PARAMS:
      u = number of the upper energy level (0=ground state)
      d = number of the lower energy level
      x = array of occupation numbers
      g = weights of all levels
      E = energy of all levels (eV) '''
def calc_T_ex(u, d, x, g, E):
    if x[u]==0.:
        return 0.0
    elif x[d]==0.:
        return 0.0 #should be inf but this is annoying to plot
    else:
        return -(E[u]-E[d])/(K_b_ev * np.log(x[u]*g[d]/(x[d]*g[u])))
    

