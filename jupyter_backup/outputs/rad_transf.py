
# coding: utf-8

# In[14]:


import numpy as np
import matplotlib.pyplot as plt
import os
import copy
from matplotlib.ticker import FormatStrFormatter as ticks
from scipy.special import erf


# In[16]:


# CONSTANTS
h_ev = 4.135667662e-15 #eV*s
h_si = 6.626070040e-34 #J*s

c_si = 2.99792458e8 #m/s
c_cgs = 2.9979245800e10 #cm/s

kb_ev = 8.6173303e-5 #eV/K                 # Boltzman constant
kb_cgs = 1.38064852e-16 # cm2 g s-2 K-1   
kb_si = 1.38064852e-23 # m2 kg s-2 K-1 or J/K  

G_si = 6.67408e-11 #m3 kg-1 s-2      # gravitational constant
G_cgs = 6.67408e-8 #cm3 g-1 s-2        

m_H_g = 1.6737236e-24    #g            # hydrogen mass
m_H_kg = 1.6737236e-27   #kg
m_p = 1.672621898e-24  #g              #mass of proton

PC = 3.0857e16 #m                       # 1 parsec
AU_cgs = 1.49597871e13 #cm              # 1 astronomical unit

year = 3.1556926e7 #s                   # 1 year

M_sun = 1.989e30 #kg                    # solar mass
L_sun_ev= 2.3892497766e+45 #eV/s        #one nominal solar luminosity
L_sun_watts = 3.826e26 #Watts = J/s     #one nominal solar luminosity

eV = 1.60217733e-19 #J
Joule = 6.241509e18 #eV


# In[17]:


T = 10.0 #K fixed for now
mu = 2.37
c_s = np.sqrt(kb_cgs*T/(mu*m_H_g))


# In[15]:


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
        E.append(float(words[1]) * c_cgs*h_ev)  # cm^-1 -> eV
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
                B[i][j] = A[i][j] * (c_cgs**2) /                     (2*h_ev * (freq[i][j])**3)  # cm2/(eV*s)
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
                        C[t][j][i] = C[t][i][j] *                             np.exp(-(E[i]-E[j])/(kb_ev*temps[t]))*g[i]/g[j]
        # add collision partner data to global array
        C_all.append(C)

    f.close()
    C_all = np.array(C_all)
    temps_all = np.array(temps_all)
    return mu, num_lvls, np.array(E), np.array(g), np.array(freq), np.array(A), np.array(B), C_all, num_partners, temps_all

