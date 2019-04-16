
# coding: utf-8

# In[1]:


get_ipython().run_cell_magic(u'time', u'', u'#CORRECT\n#what we want\n#just need to speed it up (the self-shield. iterations function is the bottleneck), and\n# then apply to all cells in the simulation\nimport numpy as np\nimport matplotlib.pyplot as plt\nimport pylab as py\nimport matplotlib.colors as cm\nimport pynbody\nfrom matplotlib.lines import Line2D\nfrom scipy.special import erf\nfrom michaels_functions import center_and_r_vir, remove_bulk_velocity\nfrom matplotlib.colors import LogNorm\nfrom matplotlib.pyplot import figure')


# In[13]:


get_ipython().run_cell_magic(u'time', u'', u'def make_pdf(s, s_bar, sigma_s):\n    pdf = (1/np.sqrt(2*np.pi*(sigma_s**2))) * (np.exp(-0.5*(((s - s_bar)/sigma_s)**2)))\n    return np.array(pdf)\n\ndef calc_integral2(s, pdf, n_H_mean, X_H2, ds):\n    integ = 0.0\n    for i in range(0, 100):\n        integ += (np.exp(s[i])*pdf[i]*X_H2[i]*ds)\n    return integ\n\ndef calc_integral(s_bar, sigma_s, smin, smax):\n    I = 0.5*np.exp(s_bar + ((sigma_s**2)/2))*( erf( (-smin+s_bar+(sigma_s**2))/(np.sqrt(2)*sigma_s) ) - erf( (-smax+s_bar+(sigma_s**2))/(np.sqrt(2)*sigma_s) ) )\n    return I\n\ndef calc_lambda_jeans(n_H, temp, m_p, K_b, G):\n    lambda_jeans = ((np.sqrt(K_b * temp / m_p)) / np.sqrt(4* np.pi * G * n_H * m_p))\n    return np.array(lambda_jeans)\n\ndef calc_n_LW(n_H, G_o, lambda_jeans, m_p):\n    kappa = 1000 * m_p\n    rad_field_outside = G_o #in solar units\n    exp_tau = np.exp(-kappa * n_H * lambda_jeans)\n    n_LW = rad_field_outside * exp_tau\n    return n_LW\n\ndef calc_n_LW_ss(n_H, n_H2, G_o, lambda_jeans, m_p):\n    kappa = 1000 * m_p\n    rad_field_outside = G_o #in solar units\n    exp_tau = np.exp(-kappa * n_H * lambda_jeans)\n    N_H2 = n_H2*lambda_jeans\n    term1 = (0.965/((1+(N_H2/5e14))**2))\n    term2 = ( (0.035/np.sqrt(1+(N_H2/5e14))) * np.exp(-1*np.sqrt(1+(N_H2/5e14))/1180) )\n    S_H2 = term1 + term2\n    n_LW_ss = rad_field_outside * exp_tau * S_H2\n    return n_LW_ss, S_H2, N_H2\n\ndef calc_X_H2(n_H, Z, n_LW):\n    DC = 1.7e-11\n    CC = 2.5e-17            #cm3 s-1\n    numerator = DC * n_LW\n    denominator = CC * Z * n_H\n    X_H2 = 1 / (2 + (numerator/denominator) )\n    return X_H2\n\ndef self_shielding_iterations(n_H, G_o, lambda_jeans, Z, m_p):\n    n_LW = calc_n_LW(n_H, G_o, lambda_jeans, m_p)\n    X_H2_a = calc_X_H2(n_H, Z, n_LW)\n    n_H2_a = n_H * X_H2_a\n\n    n_LW_1, S_H2_1, N_H2_1 = calc_n_LW_ss(n_H, n_H2_a, G_o, lambda_jeans, m_p)\n    X_H2_1 = calc_X_H2(n_H, Z, n_LW_1)\n    n_H2_1 = n_H * X_H2_1\n\n    n_LW_2, S_H2_2, N_H2_2 = calc_n_LW_ss(n_H, n_H2_1, G_o, lambda_jeans, m_p)\n    X_H2_2 = calc_X_H2(n_H, Z, n_LW_2)\n    n_H2_2 = n_H * X_H2_2\n\n    n_LW_3, S_H2_3, N_H2_3 = calc_n_LW_ss(n_H, n_H2_2, G_o, lambda_jeans, m_p)\n    X_H2_3 = calc_X_H2(n_H, Z, n_LW_3)\n    n_H2_3 = n_H * X_H2_3\n\n    n_LW_4, S_H2_4, N_H2_4 = calc_n_LW_ss(n_H, n_H2_3, G_o, lambda_jeans, m_p)\n    X_H2_4 = calc_X_H2(n_H, Z, n_LW_4)\n    n_H2_4 = n_H * X_H2_4\n\n    n_LW_5, S_H2_5, N_H2_5 = calc_n_LW_ss(n_H, n_H2_4, G_o, lambda_jeans, m_p)\n    X_H2_5 = calc_X_H2(n_H, Z, n_LW_5)\n    n_H2_5 = n_H * X_H2_5\n\n    n_LW_6, S_H2_6, N_H2_6 = calc_n_LW_ss(n_H, n_H2_5, G_o, lambda_jeans, m_p)\n    X_H2_6 = calc_X_H2(n_H, Z, n_LW_6)\n    n_H2_6 = n_H * X_H2_6\n\n    n_LW_7, S_H2_7, N_H2_7 = calc_n_LW_ss(n_H, n_H2_6, G_o, lambda_jeans, m_p)\n    X_H2_7 = calc_X_H2(n_H, Z, n_LW_7)\n    n_H2_7 = n_H * X_H2_7\n\n    n_LW_8, S_H2_8, N_H2_8 = calc_n_LW_ss(n_H, n_H2_7, G_o, lambda_jeans, m_p)\n    X_H2_8 = calc_X_H2(n_H, Z, n_LW_8)\n    n_H2_8 = n_H * X_H2_8\n\n    n_LW_9, S_H2_9, N_H2_9 = calc_n_LW_ss(n_H, n_H2_8, G_o, lambda_jeans, m_p)\n    X_H2_9 = calc_X_H2(n_H, Z, n_LW_9)\n    n_H2_9 = n_H * X_H2_9\n\n    n_LW_10, S_H2_10, N_H2_10 = calc_n_LW_ss(n_H, n_H2_9, G_o, lambda_jeans, m_p)\n    X_H2_10 = calc_X_H2(n_H, Z, n_LW_10)\n    n_H2_10 = n_H * X_H2_10\n\n    n_LW_ss, S_H2, N_H2 = calc_n_LW_ss(n_H, n_H2_10, G_o, lambda_jeans, m_p)\n    X_H2 = calc_X_H2(n_H, Z, n_LW_ss)\n    n_H2 = n_H * X_H2\n\n    return n_LW, X_H2_a, n_H2_a, n_LW_ss, S_H2, N_H2, X_H2, n_H2\n')


# In[17]:


get_ipython().run_cell_magic(u'time', u'', u'def inside_loop(M, n_H_mean, m_p, K_b, G, temp, Z_z, G_o, s, pdf, lambda_jeans, X_H2, n_H, n_LW, n_LW_ss, S_H2, N_H2, X_H2_a, n_H2_a, n_H2):\n    sigma_s = np.sqrt(np.log(1 + ((0.3 * M)**2)))\n    s_bar = -0.5*(sigma_s**2)\n    smin = -7*sigma_s + s_bar\n    smax = 7*sigma_s + s_bar\n    ds = (smax - smin)/100\n    for i in range(0, 100):\n        s[i] = smin + i*ds\n    pdf = make_pdf(s, s_bar, sigma_s)\n    n_H = n_H_mean * np.exp(s)\n    lambda_jeans = calc_lambda_jeans(n_H, temp, m_p, K_b, G)\n    for i in range(0, 100):\n        n_LW[i], X_H2_a[i], n_H2_a[i], n_LW_ss[i], S_H2[i], N_H2[i], X_H2[i], n_H2[i] = self_shielding_iterations(n_H[i], G_o, lambda_jeans[i], Z_z, m_p)\n    X_H2_bar = calc_integral2(s, pdf, n_H_mean, X_H2, ds)\n\n    return X_H2_bar')


# In[3]:


get_ipython().run_cell_magic(u'time', u'', u'path = "bulk1/data_2/hydro_59/output/"\ndata = pynbody.load(path + "output_00050")\naexp = data.properties[\'a\']\ndata.physical_units()\nr_vir = center_and_r_vir(data, aexp, path)\nremove_bulk_velocity(data)\nr_e = 0.1 * r_vir\nsph_5 = pynbody.filt.Sphere(radius = \'5 kpc\') #%(r_e*1.4))\nregion = data[sph_5]\nrho = region.gas["rho"].in_units("m_p cm^-3")\nf = open(data.filename + "/info_"+data.filename[-5:]+".txt","r")\nlines = f.readlines()\nf.close()\nfor line in lines:\n    if line[0:13]=="unit_l      =":\n        print line[:-1]\n        unit_l = float(line[14:-1])\n    if line[0:13]=="unit_d      =":\n        print line[:-1]\n        unit_d = float(line[14:-1])\n    if line[0:13]=="unit_t      =":\n        print line[:-1]\n        unit_t = float(line[14:-1])\n    if line[0:13]=="omega_b     =":\n        print line[:-1]\n        omega_b = float(line[14:-1])\nturb = np.sqrt( region.g["turb"] * 2./3. ) * unit_l / unit_t / 1e5\nturb = pynbody.array.SimArray(turb*1e5, units = "cm s**-1")\nc_s = np.sqrt(region.gas["p"] / region.gas["rho"])\nc_s = c_s.in_units(\'cm s**-1\')\nmach_no = turb / c_s\nregion.g["mach"] = mach_no.in_units("1")\nm_p_1 = pynbody.array.SimArray(1.0, pynbody.units.m_p)\nn_H_mean_sim = rho / m_p_1')


# In[18]:


get_ipython().run_cell_magic(u'time', u'', u'#mach_no_arr = mach_no\nmach_no_arr = np.logspace(-2, 2, 40)\n#n_H_mean_arr = n_H_mean_sim\nn_H_mean_arr = np.logspace(-1.5, 4, 40)\nlabel = "M "\nmin_n_H = np.log10(np.min(n_H_mean_arr))\nmax_n_H = np.log10(np.max(n_H_mean_arr))\nmin_M = np.min(mach_no_arr)\nmax_M = np.max(mach_no_arr)\nm_p = pynbody.array.SimArray(1.672621777e-24, "g")\nK_b = pynbody.array.SimArray(1.38064852e-16, "cm**2 g s**-2 K**-1")\nG = pynbody.array.SimArray(6.67259e-8, "cm**3 g**-1 s**-2")\n#temp = region.gas[\'temp\']\ntemp = pynbody.array.SimArray(10., "K")\nZ_z = 1\nG_o = 1\ns = np.zeros(100)\npdf = np.zeros(100)\nlambda_jeans = np.zeros(100)\nX_H2 = np.zeros(100)\nn_H = np.zeros(100)\nn_LW = np.zeros(100)\nn_LW_ss = np.zeros(100)\nS_H2 = np.zeros(100)\nN_H2 = np.zeros(100)\nX_H2_a = np.zeros(100)\nn_H2_a = np.zeros(100)\nn_H2 = np.zeros(100)\nX_H2_bar = np.zeros(len(n_H_mean_arr))\nfor m in range(0, len(mach_no_arr)):\n    M = mach_no_arr[m]\n    n_H_mean = n_H_mean_arr[m]\n    #X_H2_bar[m] = \n    inside_loop(M, n_H_mean, m_p, K_b, G, temp, Z_z, G_o, s, pdf, lambda_jeans, X_H2, n_H, n_LW, n_LW_ss, S_H2, N_H2, X_H2_a, n_H2_a, n_H2 )\nnp.save(\'X_H2_bar_small/unt2_X_H2_bar.npy\', X_H2_bar)\nnp.save(\'X_H2_bar_small/unt2_mach_no_arr.npy\', mach_no_arr)\nnp.save(\'X_H2_bar_small/unt2_n_H_mean_arr.npy\', n_H_mean_arr)')


# In[ ]:


region.gas["X_H2_bar"] = X_H2_bar


# In[ ]:


region.gas["n_H_mean_arr"] = n_H_mean_arr


# In[10]:


min_X = np.min(X_H2_bar)
max_X = np.max(X_H2_bar)


# In[12]:


get_ipython().run_cell_magic(u'time', u'', u'plt.figure(figsize=(9,5))\nhistX_H2_M_mass, yedges, xedges = np.histogram2d(X_H2_bar, np.log10(n_H_mean_arr),\n                           weights=mach_no_arr, bins=50 , range=[[min_X,max_X],[min_n_H,max_n_H]])\nhistX_H2_mass, yedges, xedges = np.histogram2d(X_H2_bar, np.log10(n_H_mean_arr),\n                           weights=None, bins=50 , range=[[min_X,max_X],[min_n_H,max_n_H]])\n\nyX_H2_M, xX_H2_M = yedges, xedges\nplt.pcolormesh(xedges, yedges, histX_H2_M_mass/histX_H2_mass, norm=LogNorm(), vmin=min_M, vmax=max_M)\nplt.colorbar(label=r"$\\mathcal{M}$")\nplt.xlabel(\'log(n_H_mean)\')\nplt.ylabel(\'X_H2_bar\')\nplt.grid(b=True, which=\'both\', axis=\'both\')\nplt.title(\'log(n_H_mean) vs X_H2_bar - M=varied, Z=1, G_o=1\')\nplt.savefig(\'X_H2_bar_small/unt2_log(n_H_mean)vsX_H2_bar--M.png\')\n')

