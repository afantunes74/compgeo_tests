#!/usr/bin/env python
# coding: utf-8

# (ch04_nb04)=
# # Angles, intersections, and poles

# This notebook illustrates the use of the functions in the module `angles` to solve several interesting problems. Let's start with the following problem: Two limbs of a chevron fold (A and B) have orientations (strike/dip, RHR) as follows:
# 
# Limb A = 120/40
# 
# Limb B = 250/60
# 
# Determine: (a) the trend and plunge of the hinge line of the fold, (b) the rake of the hinge line in limb A, (c) the rake of the hinge line in limb B. 

# In[1]:


import numpy as np
pi = np.pi

# Import functions 
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.angles import angle_bw_lines
from compgeo.angles import int_bw_planes


# In[2]:


# Strike and dip of the limbs in radians
str1, dip1 = np.radians([120, 40])
str2, dip2 = np.radians([250, 60])

# (a) Chevron folds have planar limbs. The hinge
# of the fold is the intersection of the limbs
htrd, hplg = int_bw_planes(str1,dip1,str2,dip2)
print("Hinge trend = {:.1f}, plunge {:.1f}"
      .format(htrd*180/pi,hplg*180/pi))

# The rake of the hinge on either limb is the angle 
# between the hinge and the strike line on the limb. 
# This line is horizontal and has plunge = 0
plg = 0

# (b) For the SW dipping limb
ang = angle_bw_lines(str1,plg,htrd,hplg)
print("Rake of hinge in SW dipping limb = {:.1f} E"
      .format(ang*180/pi))

# (c) And for the NW dipping limb
ang = angle_bw_lines(str2,plg,htrd,hplg)
print("Rake of hinge in NW dipping limb = {:.1f} W"
      .format(ang*180/pi))


# Let's do another problem: A quarry has two walls, one trending 002 and the other 135. The apparent dip of bedding on the faces are 40N and 30 SE respectively. Calculate the strike and dip of bedding.

# In[3]:


# Import function
from compgeo.angles import plane_from_app_dips

# The apparent dips are just two lines on bedding
# These lines have orientations:
trd1, plg1 = np.radians([2, 40])
trd2, plg2 = np.radians([135, 30])

# Calculate bedding from these two apparent dips
strike, dip = plane_from_app_dips(trd1,plg1,trd2,plg2)
print("Bedding strike = {:.1f}, dip {:.1f}"
      .format(strike*180/pi,dip*180/pi))


# And the final problem: Slickenside lineations trending 074 occur on a fault with orientation 300/50 (RHR). Determine the plunge of these lineations and their rake in the plane of the fault.

# In[4]:


# The lineation on the fault is just the intersection
# of a vertical plane with a strike equal to
# the trend of the lineation, and the fault
str1, dip1 = np.radians([74, 90])
str2, dip2 = np.radians([300, 50])

# Find the intersection of these two planes which is
# the lineation on the fault
ltrd, lplg = int_bw_planes(str1,dip1,str2,dip2)
print("Slickensides trend = {:.1f}, plunge {:.1f}"
      .format(ltrd*180/pi,lplg*180/pi))

# And the rake of this lineation is the angle
# between the lineation and the strike line on the fault
plg = 0
ang = angle_bw_lines(str2,plg,ltrd,lplg)
print("Rake of slickensides = {:.1f} W".format(ang*180/pi))


# There are many interesting problems you can solve using the functions in the module `angles`. You will find more problems in the Exercises section.  

# In[ ]:




