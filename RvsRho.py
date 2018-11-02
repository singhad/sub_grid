import numpy as np
import matplotlib.pyplot as plt

#function to generate the PDF
def make_pdf(x, x_mean, s):
    pdf = (1/np.sqrt(2*np.pi*(sigma**2))) * (np.exp(-0.5*(((x - x_mean)/sigma)**2)))
    return pdf

#function to calculate Jeans length
def jeans_R(T_mean, n):
    K_b = 1.38064852e-16    # ergs K-1
    m_p = 1.672621777e-24   # g
    G = 6.67408e-8          # dyne cm^2 g^-2
    c_s = np.sqrt(K_b * T_mean / m_p)       # cm s^-1
    lambda_jeans_cm = c_s / np.sqrt(4* np.pi * G * n * m_p)        # cm
    #r_jeans_cm = lambda_jeans * 0.5    # cm    #Jeans radius = 0.5*Jeans Length
    return lambda_jeans_cm

def jeans_R_pc(T_mean, n):
    L = 3.08567758e18   # 1 parsec = 3.08567758 * 10^(18) cm
    K_b = 1.38064852e-16    # ergs K-1
    m_p = 1.672621777e-24   # g
    G = 6.67408e-8          # dyne cm^2 g^-2
    c_s = np.sqrt(K_b * T_mean / m_p)       # cm s^-1
    lambda_jeans_cm = c_s / np.sqrt(4* np.pi * G * n * m_p)
    lambda_jeans_pc = lambda_jeans_cm/L     #in parsecs
    return lambda_jeans_pc

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
    lambda_jeans_pc = []    # pc
    for i in range(1, n_range):
        n = (np.logspace(-4, 5, 100) * i) # [H] cm^-3
        x.append(np.log(n/n_mean))
        pdf.append(make_pdf(np.log(n/n_mean), x_mean, sigma))
        lambda_jeans_cm.append(jeans_R(T_mean, n))
        lambda_jeans_pc.append(jeans_R_pc(T_mean, n))

    lambda_min = np.min(lambda_jeans_cm)
    lambda_max = np.max(lambda_jeans_cm)
    print(lambda_min)
    print(lambda_max)
    #plotting log(n) vs PDF
    plt.plot(x[1], pdf[1], lw=1, color='b')
    plt.xlabel('log(n) [H/cc]')
    plt.ylabel('PDF')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n) vs PDF')
    plt.savefig('XvsPDF.png'.format(i=i))
    plt.clf()
    #plotting log(n) vs log(PDF)
    plt.plot(x[1], np.log(pdf[1]), lw=1, color='b')
    plt.xlabel('log(n) [H/cc]')
    plt.ylabel('log(PDF)')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n) vs log(PDF)')
    plt.savefig('XvslogPDF.png'.format(i=i))
    plt.clf()

    #plotting Jeans length in cm lambda_jeans_cm vs log(n)
    plt.plot(lambda_jeans_cm[1], x[1], lw=1, color='b')
    plt.ylabel('log(n) [H/cc]')
    plt.xlabel('lambda_Jeans [cm]')
    plt.grid(b=True, which='both', axis='both')
    plt.title('lambda Jeans vs log(n)')
    plt.savefig('lambda_JeansvsX.png'.format(i=i))
    plt.clf()
    #plotting log of Jeans length in cm log(lambda_jeans_cm) vs log(n)
    plt.plot(np.log(lambda_jeans_cm[1]), x[1], lw=1, color='b')
    plt.ylabel('log(n) [H/cc]')
    plt.xlabel('log(lambda_Jeans) [cm]')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(lambda_Jeans) vs log(n)')
    plt.savefig('log(lambda_Jeans)vsX.png'.format(i=i))
    plt.clf()

    #plotting log of Jeans length in pc lambda_jeans_pc vs log(n)
    plt.plot(lambda_jeans_pc[1], x[1], lw=1, color='b')
    plt.ylabel('log(n) [H/cc]')
    plt.xlabel('lambda_Jeans [pc]')
    plt.grid(b=True, which='both', axis='both')
    plt.title('lambda_Jeans vs log(n)')
    plt.savefig('lambda_Jeans_pcvsX.png'.format(i=i))
    plt.clf()
    #plotting log of Jeans length in pc log(lambda_jeans_pc) vs log(n)
    plt.plot(np.log(lambda_jeans_pc[1]), x[1], lw=1, color='b')
    plt.ylabel('log(n) [H/cc]')
    plt.xlabel('log(lambda_Jeans) [pc]')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(lambda_Jeans) vs log(n)')
    plt.savefig('log(lambda_Jeans_pc)vsX.png'.format(i=i))
    plt.clf()
