#!/usr/bin/env python
# coding: utf-8

# (ch04_nb06)=
# # Uncertainties

# Suppose that in the first problem on page 64, the uncertainty in strike is 4$^\circ$ and in dip is 2$^\circ$. This problem can be solved as follows:

# In[1]:


# Import libraries
import numpy as np
pi = np.pi
from uncertainties import ufloat

# Import functions
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.angles_u import angle_bw_lines_u
from compgeo.angles_u import int_bw_planes_u


# In[2]:


# Strike and dip of the limbs in radians
str1, dip1 = np.radians([120, 40]) # SW dipping limb
str2, dip2 = np.radians([250, 60]) # NW dipping limb

# Errors in strike and dip in radians
ustr, udip = np.radians([4, 2])

# Create the input values with uncertainties
str1 = ufloat(str1, ustr)  # str1 = str1 +/-ustr
dip1 = ufloat(dip1, udip)  # dip1 = dip1 +/-udip
str2 = ufloat(str2, ustr)  # str2 = str2 +/-ustr
dip2 = ufloat(dip2, udip)  # dip2 = dip2 +/-udip

# (a) Chevron folds have planar limbs. The hinge
# of the fold is the intersection of the limbs
htrd, hplg = int_bw_planes_u(str1,dip1,str2,dip2)
print("Hinge trend = {:.1f}, plunge {:.1f}"
      .format(htrd*180/pi,hplg*180/pi))

# The rake of the hinge on either limb is the angle 
# between the hinge and the strike line on the limb. 
# This line is horizontal and has plunge = 0
plg = ufloat(0, udip)  # plg = 0 +/-udip

# (b) For the SW dipping limb
ang = angle_bw_lines_u(str1,plg,htrd,hplg)
print("Rake of hinge in SW dipping limb = {:.1f} E"
      .format(ang*180/pi))

# (c) And for the NW dipping limb
ang = angle_bw_lines_u(str2,plg,htrd,hplg)
print("Rake of hinge in NW dipping limb = {:.1f} W"
      .format(ang*180/pi))


# In the map of Fig. 4.6, the error in East and North coordinates is 10 m, and in elevation is 5 m. What is the strike and dip of the T-S contact? 

# In[3]:


# Import function three_points_u
from compgeo.three_points_u import three_points_u

# ENU coordinates of the three points
# with uncertainties in E-N = 10, and U = 5
p1 = np.array([ufloat(509, 10), ufloat(2041, 10), 
               ufloat(400, 5)])
p2 = np.array([ufloat(1323, 10), ufloat(2362, 10), 
               ufloat(500, 5)])
p3 = np.array([ufloat(2003, 10), ufloat(2913, 10), 
               ufloat(700, 5)])

# Compute the orientation of the plane
strike, dip = three_points_u(p1,p2,p3)
print("Plane strike = {:.1f}, dip = {:.1f}"
      .format(strike*180/pi,dip*180/pi))


# In[ ]:




