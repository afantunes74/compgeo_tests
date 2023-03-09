#!/usr/bin/env python
# coding: utf-8

# (ch04_nb02)=
# # Vector operations

# In[1]:


# Import numpy
import numpy as np


# In[2]:


# Make vectors
u = np.array([1,2,3])
v = np.array([3,2,1])
print("u =", u)
print("v =", v)
# Scalar multiplication of vector
sv = 3 * u
print("3 * u =", sv)
# Sum of vectors
vsum = u + v
print("u + v =", vsum)
#  Dot product of vectors
dotp = np.dot(u,v)
print("u . v =", dotp)
# Cross product of vectors
crossp = np.cross(u,v) 
print("u x v =", crossp)


# In[ ]:




