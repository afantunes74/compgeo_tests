#!/usr/bin/env python
# coding: utf-8

# (ch08_nb04)=
# # Pure shear

# In[1]:


# Import libraries
import numpy as np

# Import function pure_shear
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.pure_shear import pure_shear


# In[2]:


# Initial points coordinates
pts = np.zeros((16,2))
pts[:,0]=[-1,-1,-1,-1,-1,-0.5,0,0.5,1,1,1,1,1,0.5,0,-0.5]
pts[:,1]=[-1,-0.5,0,0.5,1,1,1,1,1,0.5,0,-0.5,-1,-1,-1,-1]
st1 = 2.5
ninc = 10

paths,psf,fig,ax = pure_shear(pts,st1,ninc)


# In[3]:


# Run this cell if you want to save the figure
fig.savefig("pure_shear.png", dpi=300)


# In[ ]:




