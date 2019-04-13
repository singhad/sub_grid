
# coding: utf-8

# In[7]:


import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import pynbody
import os
from michaels_functions import center_and_r_vir, config_check, halo_abundance


# In[8]:


path = "bulk1/data_2/hydro_59/output/"
data = pynbody.load(path + "output_00050")

aexp = data.properties['a']
data.physical_units()

print path
print "a =", aexp
print "z =", 1./aexp -1


# In[9]:


r_vir = center_and_r_vir(data, aexp, path)


# In[10]:


r_e = 0.1 * r_vir
print r_e


# In[11]:


sph_5 = pynbody.filt.Sphere(radius = '%f kpc' %(r_e))
h1 = data[sph_5]

