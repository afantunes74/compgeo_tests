#!/usr/bin/env python
# coding: utf-8

# (ch09_nb01)=
# # Stresses around a circular hole

# In[1]:


# Import function hoop
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.hoop import hoop


# In[2]:


# Define the geometry, 50 points along radius
# 100 points around circle
geom = [50, 100]
# Define stress: sigma1 = 50, sigma3 = 25, Pf = 0 MPa 
stress = [50, 25, 0]

# Compute and plot hoop and radial stresses
shm, srm, fig, ax = hoop(geom,stress)


# In[3]:


# Run this cell if you want to save the figure
fig.savefig("hoop.png", dpi=300)


# In[ ]:




