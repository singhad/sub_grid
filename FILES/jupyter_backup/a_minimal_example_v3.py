
# coding: utf-8

# In[5]:


import numpy as np
import matplotlib.pyplot as plt
import pynbody

from michaels_functions import center_and_r_vir, remove_bulk_velocity
from matplotlib.colors import LogNorm


# In[6]:


path = "bulk1/data_2/hydro_59/output/"
data = pynbody.load(path + "output_00050")

aexp = data.properties['a']
data.physical_units()

print path
print "a =", aexp
print "z =", 1./aexp -1


# In[7]:


r_vir = center_and_r_vir(data, aexp, path)


# In[8]:


remove_bulk_velocity(data)


# In[9]:


r_e = 0.1 * r_vir
print r_e


# In[10]:


sph_5 = pynbody.filt.Sphere(radius = '%f kpc' %(r_e*1.0))
region = data[sph_5]


# In[11]:


region.gas["rho"].in_units("m_p cm^-3")


# In[12]:


region.gas["metal"]


# In[13]:


pynbody.plot.image(region.g, qty="metal", width="%f kpc"%(2.*r_e),
                   log=True, resolution=500, cmap="magma", av_z="rho",
                   vmin=1e-4, vmax=1e-1)
plt.show()


# In[14]:


pynbody.plot.image(region.gas, units='Msol kpc^-2', width="%f kpc"%(2.*r_e),
                   log=True, resolution=500, vmin=1e5, vmax=1e9, cmap="viridis")
plt.show()


# In[15]:


region.gas["X"] = region.gas["rho"].in_units("m_p cm^-3") * region.gas["metal"]


# In[16]:


pynbody.plot.image(region.g, qty="X", width="%f kpc"%(2.*r_e),
                   log=True, resolution=500, cmap="magma", av_z="rho",
                   vmin=1e-4, vmax=1e-1)

plt.show()


# In[17]:


f = open(data.filename + "/info_"+data.filename[-5:]+".txt","r")
lines = f.readlines()
f.close()

for line in lines:
    if line[0:13]=="unit_l      =":
        print line[:-1]
        unit_l = float(line[14:-1])
    if line[0:13]=="unit_d      =":
        print line[:-1]
        unit_d = float(line[14:-1])
    if line[0:13]=="unit_t      =":
        print line[:-1]
        unit_t = float(line[14:-1])
    if line[0:13]=="omega_b     =":
        print line[:-1]
        omega_b = float(line[14:-1])


# In[18]:


turb = np.sqrt( region.g["turb"] * 2./3. ) * unit_l / unit_t / 1e5
turb = pynbody.array.SimArray(turb, units = "km s**-1")


# In[19]:


cs = np.sqrt( 4./3. * region.gas["p"] / region.gas["rho"])
cs = cs.in_units("km s**-1")

M = turb / cs


# In[20]:


region.g["mach"] = M.in_units("1")


# In[21]:


pynbody.plot.image(region.gas, qty="mach", width="%f kpc"%(2*r_e), 
                   log=True, resolution=500, cmap="magma", av_z="rho",
                   vmin=1e0, vmax=1e2, show_cbar=True)

plt.show()


# histograms in 1D and 2D:
# --

# In[22]:


plt.hist(region.g["temp"], bins=np.logspace(0.0, 4.0, 50))

plt.xlabel("$T$ [K]")
plt.ylabel("counts")
plt.xscale("log")
plt.yscale("log")
plt.show()


# In[23]:


hist_rho_T, yedges, xedges = np.histogram2d(np.log10(region.gas["temp"]),
                               np.log10(region.gas["rho"].in_units("m_p cm**-3")),
                               weights=region.gas["mach"] * region.gas["mass"], bins=50 , range=[[1,6],[-3,3]])

hist_rho_T_M, yedges, xedges = np.histogram2d(np.log10(region.gas["temp"]),
                               np.log10(region.gas["rho"].in_units("m_p cm**-3")),
                               weights=region.gas["mass"], bins=50 , range=[[1,6],[-3,3]])

y_rho_T, x_rho_T = yedges, xedges


# In[24]:


#Mach number - mass weighted
plt.pcolormesh(xedges, yedges, hist_rho_T/hist_rho_T_M, norm=LogNorm(), vmin=1.0, vmax=1e3)

plt.xlabel(r"$\log n$")
plt.ylabel(r"$\log T$")
plt.colorbar(label=r"$\mathcal{M}$")

plt.show()


# In[25]:


#just the mass
plt.pcolormesh(xedges, yedges, hist_rho_T_M, norm=LogNorm(), vmin=1e4, vmax=1e9)

plt.xlabel(r"$\log n$")
plt.ylabel(r"$\log T$")
plt.colorbar(label=r"$M \, [\mathrm{M_\odot}]$")

plt.show()


# rotate face-on and side-on:
# --

# In[26]:


with pynbody.analysis.angmom.faceon(region):
    
    pynbody.plot.image(region.gas, units='Msol kpc^-2', width="%f kpc"%(2.*r_e),
                   log=True, resolution=500, vmin=1e5, vmax=1e9, cmap="viridis")
    
    plt.show()


# In[27]:


with pynbody.analysis.angmom.sideon(region):
    
    pynbody.plot.image(region.gas, units='Msol kpc^-2', width="%f kpc"%(2.*r_e),
                   log=True, resolution=500, vmin=1e5, vmax=1e9, cmap="viridis")
    
    plt.show()

