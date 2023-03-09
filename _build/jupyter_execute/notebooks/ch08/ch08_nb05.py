#!/usr/bin/env python
# coding: utf-8

# (ch08_nb05)=
# # Simple shear

# In[1]:


# Import libraries
import numpy as np

# Import function simple_shear
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.simple_shear import simple_shear


# In[2]:


# Initial points coordinates
pts = np.zeros((16,2))
pts[:,0]=[-1,-1,-1,-1,-1,-0.5,0,0.5,1,1,1,1,1,0.5,0,-0.5]
pts[:,1]=[-1,-0.5,0,0.5,1,1,1,1,1,0.5,0,-0.5,-1,-1,-1,-1]
gamma = 2.5
ninc = 10
paths,psf,fig,ax = simple_shear(pts,gamma,ninc)


# In[3]:


# Run this cell if you want to save the figure
fig.savefig("simple_shear.png", dpi=300)


# In[ ]:




