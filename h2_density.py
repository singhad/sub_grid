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
    return lambda_jeans_cm

#function to return optical depth
def calc_optical_depth(n_H, lambda_jeans_cm):
    kappa_dust = 1000       # cm^2 g^-1
    m_p = 1.672621777e-24   # g
    kappa = kappa_dust * m_p
    tau = kappa * n_H * lambda_jeans_cm
    return tau

def calc_num_LW(tau):
    rad_field_outside = 1.0     #in solar units
    num_LW = rad_field_outside * np.exp(-tau)
    return num_LW

def calc_X_H2(n_H, n_LW):
    destruction_coeff = 1.7e-11
    creation_coeff = 2.5e-17    #cm3 s-1
    Z = 1.0     #solar
    factor = (destruction_coeff * n_LW) / (creation_coeff * Z * n_H)
    X_H2 = 1 / (2 + factor)
    return X_H2

def calc_num_H2(n_H, X_H2):
    n_H2 = n_H * X_H2   #H2/cc
    return n_H2

if __name__=='__main__':
    fig, ax = plt.subplots()
    m_p = 1.672621777e-24  # g
    T_mean = 10.    #K
    mach_no = 5.
    n_H_range = 100
    n_H_mean = 100    # [H] cm^-3
    vel_disp = np.sqrt(np.log(1 + ((0.3 * mach_no)**2))) #vel_disp in pdf
    x_mean = 1
    pdf = np.zeros(100)
    x = np.zeros(100)
    lambda_jeans_cm = np.zeros(100)    # cm
    tau = np.zeros(100)     #optical depth
    n_LW = np.zeros(100)    #number of Lyman-Werner photons
    n_H2 = np.zeros(100)
    X_H2 = np.zeros(100)

    for i in range(1, n_H_range):
        n_H = (np.logspace(-4, 5, 100) * i) # [H] cm^-3
        x = np.log(n_H/n_H_mean)
        pdf = make_pdf(np.log(n_H/n_H_mean), x_mean, vel_disp)
        lambda_jeans_cm = calc_lambda_jeans(T_mean, n_H)
        tau = calc_optical_depth(n_H, lambda_jeans_cm)
        n_LW = calc_num_LW(tau)
        X_H2 = calc_X_H2(n_H, n_LW)
        n_H2 = calc_num_H2(n_H, X_H2)

    #plotting log(n) vs log(PDF)
    plt.plot(x, np.log(pdf), lw=1, color='b')
    plt.xlabel('log(n) [H/cc]')
    plt.ylabel('log(PDF)')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n) vs log(PDF)')
    plt.savefig('XvslogPDF_fromcode.png'.format(i=i))
    plt.clf()
    #plotting log(n) vs PDF
    plt.plot(x, pdf, lw=1, color='b')
    plt.xlabel('log(n) [H/cc]')
    plt.ylabel('PDF')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n) vs PDF')
    plt.savefig('XvsPDF_fromcode.png'.format(i=i))
    plt.clf()
    #plotting log of Jeans length in cm log(lambda_jeans_cm) vs log(n)
    plt.plot(np.log(lambda_jeans_cm), x, lw=1, color='b')
    plt.ylabel('log(n) [H/cc]')
    plt.xlabel('log(lambda_Jeans) [cm]')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(lambda_Jeans) vs log(n)')
    plt.savefig('log(lambda_Jeans)vsX_fromcode.png'.format(i=i))
    plt.clf()
    #plotting X_H2 vs log(n_H)
    plt.plot(x, X_H2, lw=1, color='b')
    plt.xlabel('log(n_H) [H/cc]')
    plt.ylabel('X_H2')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs X_H2')
    plt.savefig('log(n_H)vsX_H2.png'.format(i=i))
    plt.clf()
    #plotting X_H2 vs log(n_H2)
    plt.plot(np.log(n_H2), X_H2, lw=1, color='b')
    plt.xlabel('log(n_H2) [H2/cc]')
    plt.ylabel('X_H2')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H2) vs X_H2')
    plt.savefig('log(n_H2)vsX_H2.png'.format(i=i))
    plt.clf()
