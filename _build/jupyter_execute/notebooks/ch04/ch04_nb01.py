#!/usr/bin/env python
# coding: utf-8

# (ch04_nb01)=
# # Vector components, magnitude, and unit vectors

# In[1]:


# Import numpy
import numpy as np


# In[2]:


# Make vector
v = np.array([1,2,3])
print("Vector:", v)
# Magnitude of the vector
length = np.linalg.norm(v) 
print("Magnitude of the vector:", length)
# Unit vector
v_hat = v / length
print("Unit Vector:", v_hat)
# Magnitude of unit vector
length = np.linalg.norm(v_hat) 
print("Magnitude of the unit vector:", length)


# In[ ]:




