import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import lognorm

#function to generate the PDF
def make_pdf(x, x_bar, s):
    pdf = (1/(s * np.sqrt(2*np.pi))) * np.exp(-1/2*(((x - x_bar)/s)**2))
    return pdf

#function to calculate Jeans length
def jeans_R(T_mean, rho):
    K_b = 1.38064852e-23    # m2 kg s-2 K-1
    m_p = 1.672621777e-27   # kg
    G = 1.672621777e-27     # m^3 kg^-1 s^-2
    c_s = np.sqrt(K_b * T_mean / m_p)
    r_jeans = c_s / np.sqrt(G * rho)
    return r_jeans

if __name__=='__main__':
    fig, ax = plt.subplots()
    rho_bar = 2.
    sigma_L = 1.
    L = 1.
    T_mean = 5.
    rho_range = 1000
    x_bar = sigma_L**2 / 2
    s = sigma_L
    pdf = []
    rho = []
    x = []
    r_jeans= []
    for i in range(1, rho_range):
        rho = np.logspace(-4, 5, 1000) * i
        x.append(np.log(rho/rho_bar))
        pdf.append(make_pdf(np.log(rho/rho_bar), x_bar, s))
        r_jeans.append(jeans_R(T_mean, rho))

    #plotting log(rho) vs PDF
    plt.plot(x[1], pdf[1], lw=1, color='b')
    plt.xlabel('log(rho)')
    plt.ylabel('PDF')
    plt.grid(b=True, which='both', axis='both')
    plt.savefig('XvsPDF.png'.format(i=i), dpi= 'figure', bbox_inches= 'tight')
    plt.clf()
    #plotting log(rho) vs log(PDF)
    plt.semilogy(x[1], pdf[1], lw=1, color='b')
    plt.xlabel('log(rho)')
    plt.ylabel('log(PDF)')
    plt.grid(b=True, which='both', axis='both')
    plt.savefig('XvslogPDF.png'.format(i=i), dpi= 'figure', bbox_inches= 'tight')
    plt.clf()

    #plotting Jeans length r_jeans vs log(rho)
    plt.plot(r_jeans[1], x[1], lw=1, color='b')
    plt.ylabel('log(rho)')
    plt.xlabel('R_Jeans')
    plt.grid(b=True, which='both', axis='both')
    plt.savefig('R_JeansvsX.png'.format(i=i))
    plt.clf()


    plt.plot(x[1], x[1], lw=1, color='b')
    plt.ylabel('log(rho)')
    plt.xlabel('log(rho)')
    plt.grid(b=True, which='both', axis='both')
    plt.savefig('XvsX.png'.format(i=i))
    plt.clf()

    plt.loglog(pdf[1], pdf[1], lw=1, color='b')
    plt.ylabel('log(pdf)')
    plt.xlabel('log(pdf)')
    plt.grid(b=True, which='both', axis='both')
    plt.savefig('logPDFvslogPDF.png'.format(i=i))
    plt.clf()
