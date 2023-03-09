#!/usr/bin/env python
# coding: utf-8

# (ch04_nb05)=
# # Three points problem

# In[1]:


import numpy as np
pi = np.pi

# Import function three_points
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.three_points import three_points


# In[2]:


# ENU coordinates of the three points
p1 = np.array([509, 2041, 400])
p2 = np.array([1323, 2362, 500])
p3 = np.array([2003, 2913, 700])

# Compute the orientation of the plane
strike, dip = three_points(p1,p2,p3)
print("Plane strike = {:.1f}, dip = {:.1f}"
      .format(strike*180/pi,dip*180/pi))


# In[ ]:




