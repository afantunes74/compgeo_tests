#!/usr/bin/env python
# coding: utf-8

# (ch08_nb02)=
# # 2D strain from GPS data

# In[1]:


import numpy as np
pi = np.pi

# Import function grid_strain
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.grid_strain import grid_strain


# In[2]:


# Load Zhang et al. GPS data from the Tibetan plateau
# load x, y coordinates and displacements
tibet =np.loadtxt(os.path.abspath("../../data/ch8-2/tibet.txt"))
pos = tibet[:,0:2]
disp = tibet[:,2:4] 

# Rotation from Delaunay triangulation, plot stations
par = 10 * pi/180 #Minimum internal angle of triangles
cent,eps,ome,pstrain,rotc,fig,ax = grid_strain(pos,disp,0,
                                               par,2,1)


# In[3]:


# Rotation from nearest neighbor 
# Grid spacing = 75 km, neighbors 6, 
# max. distance = 150 km, plot stations
par=[75e3,6,150e3]
cent,eps,ome,pstrain,rotc,fig,ax = grid_strain(pos,disp,1,
                                               par,2,1)


# In[4]:


# Rotation from distance Weighted 
# Grid spacing = 75 km, alpha = 150 km, plot stations
par=[75e3,150e3]
cent,eps,ome,pstrain,rotc,fig,ax = grid_strain(pos,disp,2,
                                               par,2,1)


# In[5]:


# Run this cell if you want to save the figure
fig.savefig("inf_strain_gps.png", dpi=300)


# In[ ]:




