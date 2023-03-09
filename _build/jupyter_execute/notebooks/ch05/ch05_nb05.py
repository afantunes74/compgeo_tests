#!/usr/bin/env python
# coding: utf-8

# (ch05_nb05)=
# # Rotations

# An overturned bed oriented 305/60 (RHR) has sedimentary lineations which indicate the palaeocurrent direction. These pitch at 60NW, with the current flowing up the plunge. Calculate the original trend of the paleocurrents.
# 
# Besides rotating the lineations back to their pre-tilted orientation, there is an additional challenge in this problem. We need to figure out the orientation of the current lineations from their pitch on the bed. We will do this as well using a rotation.

# In[1]:


import numpy as np
pi = np.pi

# Import functions 
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.pole import pole_from_plane
from compgeo.rotate import rotate
from compgeo.zero_twopi import zero_twopi


# In[2]:


# Strike and dip of bed in radians
strike, dip = np.radians([305, 60]) 

# Pole of bed
rtrd, rplg = pole_from_plane(strike, dip)

# To find the orientation of the lineations
# rotate the strike line clockwise about the 
# pole an amount equal to the pitch

# strike line
trd, plg = strike, 0 

# rotation = pitch in radians
rot = 60 * pi/180 

# orientation of lineations
trdr, plgr = rotate(rtrd,rplg,rot,trd,plg,"a")

# Now we need to rotate the lineations about
# the strike line to their pre-tilted orientation

# The bed is overturned, so it has been rotated 
# pass the vertical. The amount of rotation
# required to restore the bed to its pre-tilted
# orientation is 180- 60 = 120 deg, and it
# should be clockwise
rot = 120 * pi/180 # rotation in radians

# rotate lineations to their pre-tilted orientation
trdl, plgl = rotate(trd,plg,rot,trdr,plgr,"a")

# The current flows up the plunge, 
# so the trend of the paleocurrents is:
trdl = zero_twopi(trdl + pi)
print("Original trend of the paleocurrents = {:.1f}"
      .format(trdl*180/pi))


# In[ ]:




