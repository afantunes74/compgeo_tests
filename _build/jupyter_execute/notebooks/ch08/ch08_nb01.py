#!/usr/bin/env python
# coding: utf-8

# (ch08_nb01)=
# # P and T axes

# In[1]:


# Import libraries
import numpy as np
pi = np.pi

# Import function pt_axes
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.pt_axes import pt_axes


# In[2]:


# Read the faults
jujuy=np.loadtxt(os.path.abspath("../../data/ch8-1/jujuy.txt"), 
                   usecols = [0,1,2,3])
fault = jujuy[:,0:2] * pi/180
slip = jujuy[:,2:4] * pi/180
sense=np.loadtxt(os.path.abspath("../../data/ch8-1/jujuy.txt"), 
                   usecols = 4, dtype = "str")

# Compute P and T axes and plot them
# Don't plot the faults and slip vectors
P,T,senseC,fig,ax= pt_axes(fault,slip,sense,0)


# In[3]:


# Run this cell if you want to save the figure
fig.savefig("pt_axes.png", dpi=300)


# In[ ]:




