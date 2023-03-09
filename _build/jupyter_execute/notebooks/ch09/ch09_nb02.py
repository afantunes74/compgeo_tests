#!/usr/bin/env python
# coding: utf-8

# (ch09_nb02)=
# # Lithospheric flexure

# In[1]:


# Import libraries
import numpy as np
# Import function flex2d
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.flex2d import flex2d


# In[2]:


# Geometry
geom = np.zeros(2)
geom[0] = 500e3 # extent = 500 km
geom[1] = 5e3 # Interval in x =  5 km
# Elastic and flexural parameters
elas = np.zeros(4)
elas[0] = 70e9 # Young Modulus = 70 GPa
elas[1] = 0.25 # Poisson"s ratio = 0.25
elas[2] = 30e3 # Elastic thickness = 30 km
elas[3] = 3300 # Density of mantle in kg/m^3
# loads
loads=np.loadtxt(os.path.abspath("../../data/ch9-2/loads.txt"))

# Compute deflection profile
w, wp, fig, ax = flex2d(geom,elas,loads)


# In[3]:


# Run this cell if you want to save the figure
fig.savefig("flex.png", dpi=300)


# In[ ]:




