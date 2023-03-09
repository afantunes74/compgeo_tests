#!/usr/bin/env python
# coding: utf-8

# (ch08_nb03)=
# # 2D finite strain from displacement data

# In[1]:


import numpy as np

# Import function grid_fin_strain
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.grid_fin_strain import grid_fin_strain


# In[2]:


# load x, y deformed coordinates and displacements
demfault = np.loadtxt(os.path.abspath(
    "../../data/ch8-3/demfault.txt"))
pos = demfault[:,0:2]
disp = demfault[:,2:4]

# Max. shear strain from nearest neighbor, Def. config. 
# Grid spacing = 0.2 m, neighbors 6, 
# max. distance = 0.4 m, plot stations
par = [0.2,6,0.4]
cent,eps,pstrain,dilat,maxsh,fig,ax = grid_fin_strain(pos,
                                        disp,1,1,par,3,1)

# Add units to the axes
ax.set_xlabel("x [m]")
ax.set_ylabel("y [m]");


# In[3]:


# Run this cell if you want to save the figure
fig.savefig("fin_strain_dem.png", dpi=300)


# In[ ]:




