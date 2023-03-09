#!/usr/bin/env python
# coding: utf-8

# (ch05_nb01)=
# # Stratigraphic thickness

# In[1]:


import numpy as np
pi = np.pi

# Import function true_thickness
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.true_thickness import true_thickness


# In[2]:


# Strike and dip of the unit in radians
strike, dip = np.radians([84.5, 22.5]) 

# ENU coordinates of the points
p1 = np.array([1147, 3329, 400]) 
p2 = np.array([1323, 2362, 500]) 
p3 = np.array([1105, 1850, 400]) 
p4 = np.array([1768, 940, 300]) 
p5 = np.array([1842, 191, 200])

# Compute the thickness of the units
thickT = true_thickness(strike,dip,p2,p1)
thickS = true_thickness(strike,dip,p3,p2)
thickR = true_thickness(strike,dip,p4,p3)
thickQ = true_thickness(strike,dip,p5,p4)
print("Thickness of unit T = {:.1f} m".format(thickT))
print("Thickness of unit S = {:.1f} m".format(thickS))
print("Thickness of unit R = {:.1f} m".format(thickR))
print("Thickness of unit Q = {:.1f} m".format(thickQ))


# In[ ]:




