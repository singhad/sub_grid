'''import numpy as np
import matplotlib.pyplot as plt
import os
from matplotlib.lines import Line2D

def make_pdf(s, s_bar, sigma_s):
    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_bar)/sigma_s)**2))) * ((np.exp(s))**0.5)
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
        lambda_jeans = np.zeros(1000)
        X_H2 = np.zeros(1000)
        n_H = np.zeros(1000)
        for i in range(0, 1000):
            s[i] = smin + i*ds
        n_H = n_H_mean * np.exp(s)
        pdf = make_pdf(s, s_bar, sigma_s)
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

def make_pdf(s, s_bar, sigma_s):
    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_bar)/sigma_s)**2)))
    return pdf

def get_evolution():
    fig, ax = plt.subplots()
    ratio_arr = np.array([0.1, 0.4, 0.7, 1.])
    label = "(t/t_ff) "
    color_arr = ['r','g','b','y']
    custom_lines = [Line2D([0], [0], color='r', lw=4),
                    Line2D([0], [0], color='g', lw=4),
                    Line2D([0], [0], color='b', lw=4),
                    Line2D([0], [0], color='y', lw=4)]
    for r in range(0, len(ratio_arr)):
        mach_no = 5
        ratio = ratio_arr[r]
        color = str(color_arr[r])
        sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
        s_bar = -0.5*(sigma_s**2)
        s_tail = np.abs(0.2 * s_bar * (ratio**(-2)))
        smin = -4*sigma_s + s_bar
        smax = 4*sigma_s + s_bar
        ds = (smax - smin)/1000
        value = mach_no
        n_H_mean = 1e4
        s = np.zeros(1000)
        pdf = np.zeros(1000)
        lambda_jeans = np.zeros(1000)
        X_H2 = np.zeros(1000)
        n_H = np.zeros(1000)
        for i in range(0, 1000):
            s[i] = smin + i*ds
        n_H = n_H_mean * np.exp(s) * s_tail
        pdf = make_pdf(s, s_bar, sigma_s)
        plt.plot(np.log10(n_H), pdf, color=color)
    plt.xlabel('log(n_H)')
    plt.ylabel('pdf')
    plt.grid(b=True, which='both', axis='both')
    plt.title('log(n_H) vs pdf : varying (t/t_ff)')
    ax.legend(  custom_lines,
                [   label + '= 0.1',
                    label + '= 0.4',
                    label + '= 0.7',
                    label + '= 1.0'  ],
                loc = 'lower right'
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
