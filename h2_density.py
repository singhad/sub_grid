import numpy as np
import matplotlib.pyplot as plt

#function to generate the PDF
def make_pdf(x, x_mean, s):
    pdf = (1/np.sqrt(2*np.pi*(sigma**2))) * (np.exp(-0.5*(((x - x_mean)/sigma)**2)))
    return pdf

#function to calculate Jeans length
def get_lambda_jeans(T_mean, n):
    K_b = 1.38064852e-16    # ergs K-1
    m_p = 1.672621777e-24   # g
    G = 6.67408e-8          # dyne cm^2 g^-2
    c_s = np.sqrt(K_b * T_mean / m_p)       # cm s^-1
    lambda_jeans_cm = c_s / np.sqrt(4* np.pi * G * n * m_p)        # cm
    #r_jeans_cm = lambda_jeans * 0.5    # cm    #Jeans radius = 0.5*Jeans Length
    return lambda_jeans_cm

if __name__=='__main__':
    fig, ax = plt.subplots()
    m_p = 1.672621777e-24  # g
    T_mean = 10.           #K
    mach_no = 5.
    n_range = 1000000
    n_mean = 100    # [H] cm^-3

    sigma = np.sqrt(np.log(1 + ((0.3 * mach_no)**2))) #sigma in pdf
    x_mean = 1
    pdf = []
    n = []    # [H] cm^-3
    x = []
    lambda_jeans_cm = []     # cm
    for i in range(1, n_range):
        n = (np.logspace(-4, 5, 100) * i) # [H] cm^-3
        x.append(np.log(n/n_mean))
        pdf.append(make_pdf(np.log(n/n_mean), x_mean, sigma))
        lambda_jeans_cm.append(get_lambda_jeans(T_mean, n))
