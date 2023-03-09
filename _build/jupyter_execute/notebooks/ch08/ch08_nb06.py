#!/usr/bin/env python
# coding: utf-8

# (ch08_nb06)=
# # General shear

# In[1]:


# Import libraries
import numpy as np

# Import function general_shear
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.general_shear import general_shear


# In[2]:


# Initial points coordinates
pts = np.zeros((16,2))
pts[:,0]=[-1,-1,-1,-1,-1,-0.5,0,0.5,1,1,1,1,1,0.5,0,-0.5]
pts[:,1]=[-1,-0.5,0,0.5,1,1,1,1,1,0.5,0,-0.5,-1,-1,-1,-1]
st1 = 2.0
gamma = 0.5
ninc = 10

# Max. finite stretch parallel to shear direction
kk = 0
paths1,wk1,psf1,fig,ax= general_shear(pts,st1,gamma,kk,ninc)
print("Wk = {:.4f}".format(wk1))


# In[3]:


# Max. finite stretch perpendicular to shear direction
kk = 1
paths2,wk2,pfs2,fig,ax= general_shear(pts,st1,gamma,kk,ninc)
print("Wk = {:.4f}".format(wk2))


# In[4]:


# Run this cell if you want to save the figure
fig.savefig("general_shear.png", dpi=300)


# In[ ]:




