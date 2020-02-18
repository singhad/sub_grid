
# coding: utf-8

# In[16]:


import numpy as np
import matplotlib.pyplot as plt
import pynbody

from michaels_functions import center_and_r_vir, remove_bulk_velocity


# In[18]:


path = "bulk1/data_2/hydro_59/output/"
data = pynbody.load(path + "output_00050")

aexp = data.properties['a']
data.physical_units()

print path
print "a =", aexp
print "z =", 1./aexp -1


# In[19]:


r_vir = center_and_r_vir(data, aexp, path)


# In[20]:


remove_bulk_velocity(data)


# In[21]:


r_e = 0.1 * r_vir
print r_e


# In[22]:


sph_5 = pynbody.filt.Sphere(radius = '%f kpc' %(r_e*1.0))
region = data[sph_5]


# In[23]:


region.gas["rho"].in_units("m_p cm^-3")


# In[24]:


region.gas["metal"]


# In[25]:


pynbody.plot.image(region.g, qty="metal", width="%f kpc"%(2.*r_e),
                   log=True, resolution=500, cmap="magma", av_z="rho",
                   vmin=1e-4, vmax=1e-1)
plt.show()


# In[26]:


pynbody.plot.image(region.gas, units='Msol kpc^-2', width="%f kpc"%(2.*r_e),
                   log=True, resolution=500, vmin=1e5, vmax=1e9, cmap="viridis")
plt.show()


# In[27]:


region.gas["X"] = region.gas["rho"].in_units("m_p cm^-3") * region.gas["metal"]


# In[28]:


pynbody.plot.image(region.g, qty="X", width="%f kpc"%(2.*r_e),
                   log=True, resolution=500, cmap="magma", av_z="rho",
                   vmin=1e-4, vmax=1e-1)

plt.show()


# In[29]:


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


# In[30]:


turb = np.sqrt( region.g["turb"] * 2./3. ) * unit_l / unit_t / 1e5
turb = pynbody.array.SimArray(turb, units = "km s**-1")


# In[31]:


cs = np.sqrt( 4./3. * region.gas["p"] / region.gas["rho"])

M = turb / cs


# In[32]:


region.g["mach"] = M.in_units("1")


# In[33]:


pynbody.plot.image(region.gas, qty="mach", width="%f kpc"%(2*r_e), 
                   log=True, resolution=500, cmap="magma", av_z="rho",
                   vmin=1e0, vmax=1e2, show_cbar=True)

plt.show()


# In[34]:


with pynbody.analysis.angmom.faceon(region):
    
    pynbody.plot.image(region.gas, units='Msol kpc^-2', width="%f kpc"%(2.*r_e),
                   log=True, resolution=500, vmin=1e5, vmax=1e9, cmap="viridis")
    
    plt.show()


# In[35]:


with pynbody.analysis.angmom.sideon(region):
    
    pynbody.plot.image(region.gas, units='Msol kpc^-2', width="%f kpc"%(2.*r_e),
                   log=True, resolution=500, vmin=1e5, vmax=1e9, cmap="viridis")
    
    plt.show()

