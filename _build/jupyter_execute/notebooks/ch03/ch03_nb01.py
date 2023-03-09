#!/usr/bin/env python
# coding: utf-8

# (ch03_nb01)=
# # Plotting lines and poles in a stereonet

# In[1]:


# Import libraries
import numpy as np
import matplotlib.pyplot as plt

# Import functions pole_from_plane and
# st_coord_line
import os, sys
sys.path.append(os.path.abspath("../../"))

from compgeo.pole import pole_from_plane
from compgeo.st_coord_line import st_coord_line


# In[2]:


# Make a figure
fig, ax = plt.subplots()

# Plot the following four lines (trend and plunge)
# on an equal angle or equal area stereonet
# lines are in radians
lines = np.radians([[30, 30],[120, 45],[210, 65],[280, 15]])

# Plot the primitive of the stereonet
r = 1; # unit radius
th = np.radians(np.arange(0,361,1))
x = r * np.cos(th)
y = r * np.sin(th)
ax.plot(x,y,"k")
# Plot center of circle
ax.plot(0,0,"k+")
# Make axes equal and remove them
ax.axis("equal")
ax.axis("off")

# Find the coordinates of the lines in the
# equal angle or equal area stereonet
nrow, ncol = lines.shape
eq_angle = np.zeros((nrow, ncol))
eq_area = np.zeros((nrow, ncol))

for i in range(nrow):
    # Equal angle coordinates
    eq_angle[i,0], eq_angle[i,1] = st_coord_line(lines[i,0],
                                    lines[i,1],0) 
    # Equal area coordinates
    eq_area[i,0], eq_area[i,1] = st_coord_line(lines[i,0],
                                    lines[i,1],1)
    
# Plot the lines
# Equal angle as black dots
ax.plot(eq_angle[:,0],eq_angle[:,1],"ko")
# Equal area as red dots
ax.plot(eq_area[:,0],eq_area[:,1],"ro")
plt.show()


# In[3]:


# Plot the following four planes (strike and dip, RHR)
# as poles on an equal angle or equal area stereonet
# planes are in radians
planes = np.radians([[0, 30], [90, 50], 
                     [180, 15], [270, 65]])

# make a figure
fig, ax = plt.subplots()

# Plot the primitive of the stereonet
ax.plot(x,y,"k")
# Plot center of circle
ax.plot(0,0,"k+")
# Make axes equal and remove them
ax.axis("equal")
ax.axis("off")

# Find the coordinates of the poles to the planes in the
# equal angle or equal area stereonet
for i in range(nrow):
    # Compute pole of plane
    trd, plg = pole_from_plane(planes[i,0],planes[i,1])
    # Equal angle coordinates
    eq_angle[i,0], eq_angle[i,1] = st_coord_line(trd,plg,0) 
    # Equal area coordinates
    eq_area[i,0], eq_area[i,1] = st_coord_line(trd,plg,1)

# Plot the poles
# Equal angle as black asterisks
ax.plot(eq_angle[:,0],eq_angle[:,1],"k*")
# Equal area as red asterisks
ax.plot(eq_area[:,0],eq_area[:,1],"r*")
plt.show()


# In[ ]:




