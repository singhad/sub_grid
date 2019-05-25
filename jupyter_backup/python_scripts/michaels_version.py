#!/usr/bin/env python

import numpy as np
import pynbody
from michaels_functions import (center_and_r_vir, remove_bulk_velocity,
                                read_unit_from_info)


def make_pdf(s, s_bar, sigma_s):
    return ((1./np.sqrt(2*np.pi*(sigma_s**2))) *
            (np.exp(-0.5*(((s - s_bar)/sigma_s)**2))))


def calc_lambda_jeans(n_H, T_mean, m_p, K_b):
    return np.asarray(np.sqrt(K_b * T_mean/m_p) /
                      np.sqrt(4.*np.pi * G * n_H * m_p))


def calc_n_LW(n_H, G_o, lambda_jeans, Z, m_p):
    kappa = 1000 * m_p * Z
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    return G_o * exp_tau


def calc_X_H2(n_H, Z, n_LW):
    DC = 1.7e-11
    CC = 2.5e-17            # cm3 s-1
    numerator = DC * n_LW
    denominator = CC * Z * n_H
    X_H2 = 1. / (2. + (numerator/denominator))
    return np.asarray(X_H2)


def calc_n_LW_ss(n_H, n_H2, Z, G_o, lambda_jeans, m_p):
    kappa = 1000 * m_p * Z
    rad_field_outside = G_o  # in solar units
    exp_tau = np.exp(-kappa * n_H * lambda_jeans)
    N_H2 = n_H2 * lambda_jeans
    term1 = 0.965/((1+(N_H2/5e14))**2)
    term2 = ((0.035/np.sqrt(1+(N_H2/5e14))) *
             np.exp(-1*np.sqrt(1+(N_H2/5e14))/1180))
    S_H2 = term1 + term2
    n_LW_ss = rad_field_outside * exp_tau * S_H2
    return n_LW_ss


def self_shielding_iterations(n_H, G_o, lambda_jeans, Z, m_p):
    ctr = 16
    n_LW = calc_n_LW(n_H, G_o, lambda_jeans, Z, m_p)
    X_H2 = calc_X_H2(n_H, Z, n_LW)
    n_H2 = n_H * X_H2
    n_H2_ss = n_H2
    for _ in range(ctr):
        n_LW_ss = calc_n_LW_ss(n_H, n_H2_ss, Z, G_o, lambda_jeans, m_p)
        X_H2_ss = calc_X_H2(n_H, Z, n_LW_ss)
        """
        if (np.sum(np.square(n_H2_ss - n_H * X_H2_ss)) < 1e-5):
            n_H2_ss = n_H * X_H2_ss
            break
        """
        n_H2_ss = n_H * X_H2_ss
    return n_LW, n_H2, n_LW_ss, X_H2_ss, n_H2_ss


def calc_integral(s, pdf, X, ds):
    return np.sum(np.exp(s) * pdf * X * ds)


def calc_X_CO(n_H, n_H2, n_LW):
    rate_CHX = 5.0e-10 * n_LW
    rate_CO = 1.0e-10 * n_LW
    x0 = 2.0e-4
    k0 = 5.0e-16  # cm3 s-1
    k1 = 5.0e-10  # cm3 s-1
    factor_beta = rate_CHX/(n_H*k1*x0)
    beta = 1./(1.+factor_beta)
    factor_CO = rate_CO/(n_H2*k0*beta)
    X_CO = 1./(1.+factor_CO)
    return np.asarray(X_CO)


def inside_loop(mach_no, n_H_mean, Z, G_o, T_mean, m_p, K_b):
    sigma_s = np.sqrt(np.log(1 + ((0.3 * mach_no)**2)))
    s_bar = -0.5*(sigma_s**2)
    smin = -7*sigma_s + s_bar
    smax = 7*sigma_s + s_bar

    s = np.linspace(smin, smax, num=100, endpoint=False)
    ds = np.diff(s)[0]

    n_H = n_H_mean * np.exp(s)
    pdf = make_pdf(s, s_bar, sigma_s)
    lambda_jeans = calc_lambda_jeans(n_H, T_mean, m_p, K_b)

    n_LW, n_H2, n_LW_ss, X_H2_ss, n_H2_ss = self_shielding_iterations(
        n_H, G_o, lambda_jeans, Z, m_p)

    X_H2_bar = 2.0 * calc_integral(s, pdf, X_H2_ss, ds)
    X_CO = calc_X_CO(n_H, n_H2, n_LW)
    X_CO_bar = calc_integral(s, pdf, X_CO, ds)

    return X_H2_bar, X_CO_bar


if __name__ == '__main__':
    run = "hydro_59"
    out = "output_00050"
    path = "../../bulk1/data_2/" + run + "/output/"
    data = pynbody.load(path + out)
    aexp = data.properties['a']
    data.physical_units()

    r_vir = center_and_r_vir(data, aexp, path)
    remove_bulk_velocity(data)
    r_e = 0.1 * r_vir

    sph_5 = pynbody.filt.Sphere(radius='%f kpc' % (r_e*1.4))
    region = data[sph_5]

    omega_b, unit_l, unit_d, unit_t = read_unit_from_info(data)

    m_p = pynbody.array.SimArray(1.672621777e-24, "g")
    K_b = pynbody.array.SimArray(1.38064852e-16, "cm**2 g s**-2 K**-1")
    G = pynbody.array.SimArray(6.67259e-8, "cm**3 g**-1 s**-2")
    T_mean = pynbody.array.SimArray(10., "K")

    rho = region.gas["rho"].in_units("m_p cm**-3")

    turb = np.sqrt(region.g["turb"] * 2./3.) * unit_l / unit_t / 1e5
    turb = pynbody.array.SimArray(turb*1e5, units="cm s**-1")

    temperature = region.g["temp"]
    c_s_arr = np.sqrt(K_b * temperature / m_p)

    mach_no_sim = np.array(turb / c_s_arr)

    m_p_1 = pynbody.array.SimArray(1.0, pynbody.units.m_p)
    n_H_mean_sim = np.array(rho / m_p_1)

    Z_arr = region.g["metal"]/0.02
    G_o = 1

    mask_relevant = np.logical_and(temperature < 1e4, n_H_mean_sim > 1e-2)
    original_order = np.arange(len(rho))
    rel_cells = original_order[mask_relevant]
    X_H2_bar_cells = np.zeros(len(rho))
    X_CO_bar_cells = np.zeros(len(rho))

    for cell in rel_cells:
        mach_no = mach_no_sim[cell]
        n_H_mean = n_H_mean_sim[cell]
        Z = Z_arr[cell]

        X_H2_bar, X_CO_bar = inside_loop(
            mach_no, n_H_mean, Z, G_o, T_mean, m_p, K_b)
        X_H2_bar_cells[cell] = X_H2_bar
        X_CO_bar_cells[cell] = X_CO_bar

    np.save('outputs/X_H2_bar_' + run + '_' + out + '.npy', X_H2_bar_cells)
    np.save('outputs/X_CO_bar_' + run + '_' + out + '.npy', X_CO_bar_cells)
    np.save('outputs/mach_no_arr_' + run + '_' + out + '.npy', mach_no_sim)
    np.save('outputs/n_H_mean_arr_' + run + '_' + out + '.npy', n_H_mean_sim)
    np.save('outputs/Z_arr_' + run + '_' + out + '.npy', Z_arr)
    np.save('outputs/T_' + run + '_' + out + '.npy', temperature)
