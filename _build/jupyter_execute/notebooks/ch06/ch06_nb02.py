#!/usr/bin/env python
# coding: utf-8

# (ch06_nb02)=
# # Best-fit plane

# In[1]:


# Import libraries
import numpy as np
pi = np.pi

# Import function fit_plane
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.fit_plane import fit_plane


# In[2]:


# Read the points on the contact
# Coordinates are UTM (ENU) in meters
jske = np.loadtxt(os.path.abspath('../../data/ch6-2/jske.txt'))

# Compute best-fit plane
strike, dip, stdev = fit_plane(jske)

# Print strike and dip of plane
print('Strike = {:.1f}, Dip = {:.1f}'
      .format(strike*180/pi,dip*180/pi))

# Print standard deviation of the distance of each point
# from the best-fit plane
print('Standard deviation = {:.1f} m'.format(stdev))


# In[ ]:




