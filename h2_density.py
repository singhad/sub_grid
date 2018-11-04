import numpy as np
import matplotlib.pyplot as plt

#function to generate the PDF
def make_pdf(x, x_mean, s):
    pdf = (1/np.sqrt(2*np.pi*(vel_disp**2))) * (np.exp(-0.5*(((x - x_mean)/vel_disp)**2)))
    return pdf

#function to calculate Jeans length
def calc_lambda_jeans(T_mean, n_H):
    K_b = 1.38064852e-16    # ergs K-1
    m_p = 1.672621777e-24   # g
    G = 6.67408e-8          # dyne cm^2 g^-2
    c_s = np.sqrt(K_b * T_mean / m_p)       # cm s^-1
    lambda_jeans_cm = c_s / np.sqrt(4* np.pi * G * n_H * m_p)        # cm
    #r_jeans_cm = lambda_jeans * 0.5    # cm    #Jeans radius = 0.5*Jeans Length
    return lambda_jeans_cm

"""#function to return optical depth
def calc_optical_depth(sigma_H2, n_H, lambda_jeans_cm):
    return sigma_H2 * n_H * lambda_jeans_cm"""

"""def calc_num_LW(tau):
    local_rad_field = 1.0  #solar units
    num_LW = local_rad_field * np.exp(-tau)
    return num_LW"""

def calc_X_H2(n_H, Z, n_LW):
    k2 = 2.5e-17        #cm3 s-1
    rate_H2 = 1.7e-11 * n_LW        #s-1
    factor = rate_H2/(2*n_H*k2*Z)
    X_H2 = 0.5 * 1./(1.+factor)
    return X_H2

def calc_num_H2(n_H, X_H2):
    return n_H * X_H2 #H2/cc

if __name__=='__main__':
    fig, ax = plt.subplots()
    m_p = 1.672621777e-24  # g
    T_mean = 10.           #K
    mach_no = 5.
    n_H_range = 1000000
    n_H_mean = 100    # [H] cm^-3
    vel_disp = np.sqrt(np.log(1 + ((0.3 * mach_no)**2))) #vel_disp in pdf
    x_mean = 1
    pdf = []
    n_H = []    # [H] cm^-3
    x = []
    lambda_jeans_cm = []     # cm
    tau = []    #optical depth
    sigma_H2 = 1000 * m_p       # assuming solar metallicity #cm^2
    n_LW = []     #number of Lyman-Werner photons
    local_rad_field = 1.0  #solar units
    n_H2 = []       #[H]/cc
    Z = 1. # solar metallicity
    X_H2 = []
    for i in range(1, n_H_range):
        n_H = (np.logspace(-4, 5, 100) * i) # [H] cm^-3
        x.append(np.log(n_H/n_H_mean))
        pdf.append(make_pdf(np.log(n_H/n_H_mean), x_mean, vel_disp))
        lambda_jeans_cm.append(calc_lambda_jeans(T_mean, n_H))
        #tau.append(calc_optical_depth(sigma_H2, n_H, lambda_jeans_cm))
        n_LW.append(np.exp(-sigma_H2 * n_H * lambda_jeans_cm)) #tau = optical depth
        X_H2.append(calc_X_H2(n_H, Z, n_LW))
        n_H2.append(calc_num_H2(n_H, X_H2))

    #plotting X_H2 vs n_H2
    plt.plot(X_H2[1], n_H2[1], lw=1, color='b')
    plt.xlabel('X_H2')
    plt.ylabel('n_H2 [H/cc]')
    plt.grid(b=True, which='both', axis='both')
    plt.title('X_H2 vs n_H2')
    plt.savefig('X_H2vsn_H2.png'.format(i=i))
    plt.clf()
