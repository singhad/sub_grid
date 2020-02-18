
# coding: utf-8

# In[2]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as cm
import pynbody
from michaels_functions import center_and_r_vir, remove_bulk_velocity
from matplotlib.colors import LogNorm
from matplotlib.pyplot import figure
import matplotlib.tri as mtri
from mpl_toolkits.mplot3d import Axes3D


# In[3]:


get_ipython().run_cell_magic(u'time', u'', u'path = "bulk1/data_2/hydro_59/output/"\ndata = pynbody.load(path + "output_00050")\naexp = data.properties[\'a\']\ndata.physical_units()\nr_vir = center_and_r_vir(data, aexp, path)\nremove_bulk_velocity(data)\nr_e = 0.1 * r_vir\nsph_5 = pynbody.filt.Sphere(radius = \'5.0 kpc\') # %(r_e*1.4))\nregion = data[sph_5]\nrho = region.gas["rho"].in_units("m_p cm**-3")\nf = open(data.filename + "/info_"+data.filename[-5:]+".txt","r")\nlines = f.readlines()\nf.close()\nfor line in lines:\n    if line[0:13]=="unit_l      =":\n        print line[:-1]\n        unit_l = float(line[14:-1])\n    if line[0:13]=="unit_d      =":\n        print line[:-1]\n        unit_d = float(line[14:-1])\n    if line[0:13]=="unit_t      =":\n        print line[:-1]\n        unit_t = float(line[14:-1])\n    if line[0:13]=="omega_b     =":\n        print line[:-1]\n        omega_b = float(line[14:-1])')


# In[4]:


get_ipython().run_cell_magic(u'time', u'', u'm_p = pynbody.array.SimArray(1.672621777e-24, "g")\nK_b = pynbody.array.SimArray(1.38064852e-16, "cm**2 g s**-2 K**-1")\nG = pynbody.array.SimArray(6.67259e-8, "cm**3 g**-1 s**-2")\nT_mean = pynbody.array.SimArray(10., "K")\nturb = np.sqrt( region.g["turb"] * 2./3. ) * unit_l / unit_t / 1e5\nturb = pynbody.array.SimArray(turb*1e5, units = "cm s**-1")\nT = region.g["temp"]\nc_s = np.sqrt(K_b * T / m_p)\nmach_no_sim = turb / c_s\nregion.g["mach"] = mach_no_sim.in_units("1")\nm_p_1 = pynbody.array.SimArray(1.0, pynbody.units.m_p)\nn_H_mean_sim = rho / m_p_1\nZ_arr = region.g["metal"]\nG_o = 1')


# In[5]:


get_ipython().run_cell_magic(u'time', u'', u"#multiplied X_H2_bar with 2 during the 1.1 iteration\nX_H2_bar = np.load('outputs/1.4/X_H2_bar_1.4.npy')\nmin_X = np.min(X_H2_bar)\nmax_X = np.max(X_H2_bar)")


# In[6]:


n_H_mean_arr = n_H_mean_sim
mach_no_arr = mach_no_sim
min_n_H = np.log10(np.min(n_H_mean_sim))
max_n_H = np.log10(np.max(n_H_mean_sim))
min_M = np.min(mach_no_sim)
max_M = np.max(mach_no_sim)
min_Z = np.min(Z_arr)
max_Z = np.max(Z_arr)


# In[7]:


get_ipython().run_cell_magic(u'time', u'', u"#multiplied X_CO_bar with 2 during the 1.1 iteration\nX_CO_bar = np.load('outputs/1.4/X_CO_bar_1.4.npy')\nmin_C = np.min(X_CO_bar)\nmax_C = np.max(X_CO_bar)")


# In[8]:


region.gas["X_H2_bar"] = X_H2_bar
region.gas["n_H_mean_arr"] = n_H_mean_arr
region.gas["X_CO_bar"] = X_CO_bar


# In[9]:


get_ipython().run_cell_magic(u'time', u'', u'plt.figure(figsize=(9,5))\nhistX_H2_M_mass, yedges, xedges = np.histogram2d(X_H2_bar, np.log10(n_H_mean_sim),\n                               weights=mach_no_sim * region.gas["mass"], bins=50 , range=[[min_X,max_X],[min_n_H,max_n_H]])\nhistX_H2_mass, yedges, xedges = np.histogram2d(X_H2_bar, np.log10(n_H_mean_sim),\n                               weights=region.gas["mass"], bins=50 , range=[[min_X,max_X],[min_n_H,max_n_H]])\n\nyX_H2_M, xX_H2_M = yedges, xedges\nplt.pcolormesh(xedges, yedges, histX_H2_M_mass/histX_H2_mass, norm=LogNorm(), vmin=min_M, vmax=max_M)\nplt.colorbar(label=r"$\\mathcal{M}$")\nplt.xlabel(\'log(n_H_mean)\')\nplt.ylabel(\'X_H2_bar\')\nplt.grid(b=True, which=\'both\', axis=\'both\')\nplt.title(\'log(n_H_mean) vs X_H2_bar - M=varied, Z=varied, G_o=1\')\nplt.savefig(\'outputs/1.4/Hist-X_H2_bar-n_H_mean-1.4.png\', dpi=300, bbox_inches=\'tight\')\nplt.show()')


# In[10]:


get_ipython().run_cell_magic(u'time', u'', u'plt.figure(figsize=(9,5))\nhistX_H2_M_mass, yedges, xedges = np.histogram2d(X_H2_bar, np.log10(n_H_mean_sim),\n                               weights=region.g["metal"] * region.gas["mass"], bins=50 , range=[[min_X,max_X],[min_n_H,max_n_H]])\nhistX_H2_mass, yedges, xedges = np.histogram2d(X_H2_bar, np.log10(n_H_mean_sim),\n                               weights=region.gas["mass"], bins=50 , range=[[min_X,max_X],[min_n_H,max_n_H]])\n\nyX_H2_M, xX_H2_M = yedges, xedges\nplt.pcolormesh(xedges, yedges, histX_H2_M_mass/histX_H2_mass, norm=LogNorm(), vmin=min_Z, vmax=max_Z)\nplt.colorbar(label="$Z$")\nplt.xlabel(\'log(n_H_mean)\')\nplt.ylabel(\'X_H2_bar\')\nplt.grid(b=True, which=\'both\', axis=\'both\')\nplt.title(\'log(n_H_mean) vs X_H2_bar - M=varied, Z=varied, G_o=1\')\nplt.savefig(\'outputs/1.4/Hist-X_H2_bar-n_H_mean_Z-1.4.png\', dpi=300, bbox_inches=\'tight\')\nplt.show()')


# In[11]:


plt.figure(figsize=(9,5))
pynbody.plot.image(region.g, qty="X_H2_bar", width='5.0 kpc',
                   log=False, resolution=1000, cmap="viridis", av_z="n_H_mean_arr",
                   vmin=5e-2, vmax=1)
plt.title("X_H2_bar-n_H_mean(cm^-3)")
plt.savefig('outputs/1.4/X_H2_bar-n_H_mean-1.4.png', dpi=300, bbox_inches='tight')
plt.show()


# In[17]:


get_ipython().run_cell_magic(u'time', u'', u'plt.figure(figsize=(9,5))\nhistX_CO_M_mass, y1edges, x1edges = np.histogram2d(X_CO_bar, np.log10(n_H_mean_sim),\n                               weights=mach_no_sim * region.gas["mass"], bins=50 , range=[[min_C,max_C],[min_n_H,max_n_H]])\nhistX_CO_mass, y1edges, x1edges = np.histogram2d(X_CO_bar, np.log10(n_H_mean_sim),\n                               weights=region.gas["mass"], bins=50 , range=[[min_C,max_C],[min_n_H,max_n_H]])\n\nyX_CO_M, xX_CO_M = y1edges, x1edges\nplt.pcolormesh(x1edges, y1edges, histX_CO_M_mass/histX_CO_mass, norm=LogNorm(), vmin=min_M, vmax=max_M)\nplt.colorbar(label="$Z$")\nplt.xlabel(\'log(n_H_mean)\')\nplt.ylabel(\'X_CO_bar\')\nplt.grid(b=True, which=\'both\', axis=\'both\')\nplt.title(\'log(n_H_mean) vs X_CO_bar - M=varied, Z=varied, G_o=1\')\nplt.savefig(\'outputs/1.4/Hist-X_CO_bar-n_H_mean-1.4.png\', dpi=300, bbox_inches=\'tight\')\nplt.show()')


# In[21]:


get_ipython().run_cell_magic(u'time', u'', u'plt.figure(figsize=(9,5))\nhistX_CO_M_mass, y1edges, x1edges = np.histogram2d(X_CO_bar, np.log10(n_H_mean_sim),\n                               weights=region.gas["metal"] * region.gas["mass"], bins=50 , range=[[min_C,max_C],[min_n_H,max_n_H]])\nhistX_CO_mass, y1edges, x1edges = np.histogram2d(X_CO_bar, np.log10(n_H_mean_sim),\n                               weights=region.gas["mass"], bins=50 , range=[[min_C,max_C],[min_n_H,max_n_H]])\n\nyX_CO_M, xX_CO_M = y1edges, x1edges\nplt.pcolormesh(x1edges, y1edges, histX_CO_M_mass/histX_CO_mass, norm=LogNorm(), vmin=min_Z, vmax=max_Z)\nplt.colorbar(label="$Z$")\nplt.xlabel(\'log(n_H_mean)\')\nplt.ylabel(\'X_CO_bar\')\nplt.grid(b=True, which=\'both\', axis=\'both\')\nplt.title(\'log(n_H_mean) vs X_CO_bar - M=varied, Z=varied, G_o=1\')\nplt.savefig(\'outputs/1.4/Hist-X_CO_bar-n_H_mean_Z-1.4.png\', dpi=300, bbox_inches=\'tight\')\nplt.show()\n')


# In[19]:


plt.figure(figsize=(9,5))
pynbody.plot.image(region.g, qty="X_CO_bar", width='5.0 kpc',
                   log=False, resolution=1000, cmap="viridis", av_z="n_H_mean_arr",
                   vmin=5e-10, vmax=1.1)
plt.title("X_CO_bar-n_H_mean(cm^-3)")
plt.savefig('outputs/1.4/X_CO_bar-n_H_mean-1.4.png', dpi=300, bbox_inches='tight')
plt.show()


# In[20]:


#insert 3D plot


# In[16]:


# %%time
# fig = plt.figure(figsize=(9,5))
# ax = fig.add_subplot(1, 2, 1, projection='3d')
# ax.plot_trisurf(X_H2_bar, np.log10(n_H_mean_sim), np.log(mach_no_sim), 
#                    norm=np.log(Z_arr), vmin=np.log(min_Z), vmax=np.log(max_Z), cmap='viridis')
# plt.colorbar(label="$Z$")
# plt.xlabel('log(n_H_mean)')
# plt.ylabel('X_H2_bar')
# plt.zlabel('Mach no.')
# plt.grid(b=True, which='both', axis='both')
# plt.title('log(n_H_mean) vs X_H2_bar - M=varied, Z=varied, G_o=1')
