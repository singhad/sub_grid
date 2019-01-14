'''import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D

def make_pdf(s, s_bar, sigma_s):
    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_bar)/sigma_s)**2))) * np.exp(s)
    return pdf

def get_evolution():
    fig, ax = plt.subplots()
    mach_no_arr = np.array([1., 5., 10., 50.])
    label = "M "
    color_arr = ['r','g','b','y']
    custom_lines = [Line2D([0], [0], color='r', lw=4),
                    Line2D([0], [0], color='g', lw=4),
                    Line2D([0], [0], color='b', lw=4),
                    Line2D([0], [0], color='y', lw=4)]
    for m in range(0, len(mach_no_arr)):
        mach_no = mach_no_arr[m]
        color = str(color_arr[m])
        sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
        s_bar = -0.5*(sigma_s**2)
        smin = -4*sigma_s + s_bar
        smax = 4*sigma_s + s_bar
        ds = (smax - smin)/1000
        value = mach_no
        n_H_mean = 1e4
        s = np.zeros(1000)
        pdf = np.zeros(1000)
        pdf2 = np.zeros(1000)
        lambda_jeans = np.zeros(1000)
        X_H2 = np.zeros(1000)
        n_H = np.zeros(1000)
        for i in range(0, 1000):
            s[i] = smin + i*ds
            n_H[i] = n_H_mean * np.exp(s[i])
            pdf[i] = make_pdf(s[i], s_bar, sigma_s)
        plt.plot(np.log10(n_H), pdf, color=color)
    plt.xlabel('log(n_H)')
    plt.ylabel('pdf')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs pdf - M=varied, Z=1, G_o=1')
    ax.legend(  custom_lines,
                [   label + '= 1',
                    label + '= 5',
                    label + '= 10',
                    label + '= 50'  ],
                loc = 'upper right'
                    )
    plt.savefig(os.path.join('log(n_H)vspdf--M.png'.format()))
    plt.clf()
    return s, smin, smax, sigma_s, n_H, pdf

if __name__=='__main__':
    path = 'for pdf_evolution'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)

    # order of variables:
    # s, smin, smax, sigma_s, n_H, pdf
    m1, m2, m3, m4, m5, m6 = get_evolution()  #varying mach_no'''


import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D
import matplotlib.gridspec as gridspec

def make_pdf(s, s_bar, sigma_s):
    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_bar)/sigma_s)**2)))
    return pdf

def get_evolution():
    fig = plt.figure()
    gs = gridspec.GridSpec(2, 4)
    ax1 = plt.subplot(gs[0, :])
    ax2 = plt.subplot(gs[1, :])
    G = 6.67408e-8          # dyne cm^2 g^-2
    m_p = 1.672621777e-24   # g
    ratio_arr = np.array([0.2, 0.4, 0.6])
    label = "ratio "
    color_arr = ['r','g','b']
    custom_lines = [Line2D([0], [0], color='k', lw=1),
                    Line2D([0], [0], color='r', lw=1, ls='--'),
                    Line2D([0], [0], color='g', lw=1, ls='--'),
                    Line2D([0], [0], color='b', lw=1, ls='--')]
    for r in range(0, len(ratio_arr)):
        mach_no = 5
        ratio = ratio_arr[r]
        color = str(color_arr[r])
        sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
        s_bar = -0.5*(sigma_s**2)
        smin = -4*sigma_s + s_bar
        smax = 4*sigma_s + s_bar
        ds = (smax - smin)/1000
        value = mach_no
        n_H_mean = 1e2
        s = np.zeros(1000)
        pdf = np.zeros(1000)
        n_H = np.zeros(1000)
        x = np.zeros(1000)
        x_prime = np.zeros(1000)
        s_prime = np.zeros(1000)
        pdf_prime = np.zeros(1000)
        n_H_prime = np.zeros(1000)
        t_ff = np.zeros(1000)
        t_ff_1 = np.zeros(1000)
        t_ff_2 = 0
        t = 0
        for i in range(0, 1000):
            s[i] = smin + i*ds
        n_H = n_H_mean * np.exp(s)
        t_ff_1 = np.sqrt((3 * np.pi)/(32 * G * n_H * m_p))
        t_ff_2 = np.sqrt((3 * np.pi)/(32 * G * n_H_mean * m_p))
        t = ratio * t_ff_2
        term0 = (1+(((t/t_ff_1)**2)*0.8614)) / ((1-((t/t_ff_1)**2))**2.8614)
        term1 = ((1+((1.8614*((t/t_ff_1)**2))/(1-((t/t_ff_1)**2))))**1)
        term2 = (-1*((t/t_ff_1)**2)/n_H) * ((2.8614/(1-((t/t_ff_1)**2))) + (0.8614/(1+(((t/t_ff_1)**2)*0.8614))))
        term3 = (1.8614/(n_H*(1-((t/t_ff_1)**2))))
        term = term3 + term2
        n_H_prime = n_H_mean * np.exp(s) * term0
        s_prime = np.log(n_H_prime/n_H_mean)
        x = np.log10(n_H/n_H_mean)
        x_prime = np.log10(n_H_prime/n_H_mean)
        pdf = make_pdf(s, s_bar, sigma_s)
        #pdf_prime = make_pdf(s_prime, s_bar, sigma_s)
        ax1.plot(x_prime, pdf, color=color, ls='--')
        ax2.plot(x_prime, np.log10(pdf), color=color, ls='--')
        ax1.set_xlim([-2, 4])
        ax2.set_xlim([-2, 4])
    ax1.plot(x, pdf, color='k', ls='-')
    ax2.plot(x, np.log10(pdf), color='k', ls='-')
    ax1.set_xlim([-2, 4])
    ax2.set_xlim([-2, 4])
    plt.xlabel('log($n_H$)')
    ax1.set_ylabel('pdf')
    ax2.set_ylabel('log(pdf)')
    ax1.grid(b=True, which='both', axis='both')
    ax2.grid(b=True, which='both', axis='both')
    fig.suptitle('log($n_H$) vs pdf and log(pdf): varying ratio')
    legend = ax1.legend(  custom_lines,
                [   label + '= 0.0',
                    label + '= 0.2',
                    label + '= 0.4',
                    label + '= 0.6'  ],
                loc = 'upper center',
                bbox_to_anchor = (0.5, 1.21),
                ncol = 4,
                fancybox = True,
                shadow = True
                    )
    plt.savefig(os.path.join('log(n_H)vspdf--ratio.png'.format()))
    plt.clf()
    return s, smin, smax, sigma_s, n_H, pdf

if __name__=='__main__':
    path = 'for pdf_evolution'
    os.makedirs(path, exist_ok=True)
    os.chdir(path)

    # order of variables:
    # s, smin, smax, sigma_s, n_H, pdf
    m1, m2, m3, m4, m5, m6 = get_evolution()  #varying mach_no
