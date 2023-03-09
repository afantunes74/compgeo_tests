#!/usr/bin/env python
# coding: utf-8

# (ch05_nb02)=
# # Stratigraphic thickness and uncertainties

# In[1]:


# Import libraries
import numpy as np
pi = np.pi
from uncertainties import ufloat

# Import function true_thickness_u
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.true_thickness_u import true_thickness_u


# In[2]:


# Strike and dip of the unit in radians
strike, dip = np.radians([84.5, 22.5]) 

# Strike and dip errors in radians
ustrike, udip = np.radians([4, 2]) 

# Create the strike and dip with uncertainties
strike = ufloat(strike, ustrike) 
dip = ufloat(dip, udip)

# ENU coordinates of the points
# with uncertainties in E-N = 10, and U = 5
p1 = np.array([ufloat(1147, 10), ufloat(3329, 10), 
               ufloat(400, 5)]) 
p2 = np.array([ufloat(1323, 10), ufloat(2362, 10), 
               ufloat(500, 5)]) 
p3 = np.array([ufloat(1105, 10), ufloat(1850, 10), 
               ufloat(400, 5)]) 
p4 = np.array([ufloat(1768, 10), ufloat(940, 10), 
               ufloat(300, 5)]) 
p5 = np.array([ufloat(1842, 10), ufloat(191, 10), 
               ufloat(200, 5)])

# Compute the thickness of the units
thickT = true_thickness_u(strike, dip, p2, p1)
thickS = true_thickness_u(strike, dip, p3, p2)
thickR = true_thickness_u(strike, dip, p4, p3)
thickQ = true_thickness_u(strike, dip, p5, p4) 
print("Thickness of unit T = {:.1f} m".format(thickT))
print("Thickness of unit S = {:.1f} m".format(thickS))
print("Thickness of unit R = {:.1f} m".format(thickR))
print("Thickness of unit Q = {:.1f} m".format(thickQ))


# In[ ]:




