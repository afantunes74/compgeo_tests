#!/usr/bin/env python
# coding: utf-8

# (ch05_nb03)=
# # Outcrop trace of a plane

# In[1]:


# Import libraries
import numpy as np
import matplotlib.pyplot as plt

# Import function outcrop_trace
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.outcrop_trace import outcrop_trace


# In[2]:


# Read the DEM grid
XG = np.loadtxt(os.path.abspath("../../data/ch5-3/XG.txt"))
YG = np.loadtxt(os.path.abspath("../../data/ch5-3/YG.txt"))
ZG = np.loadtxt(os.path.abspath("../../data/ch5-3/ZG.txt"))

# Make a figure
fig, ax = plt.subplots()

# Contour the terrain
cval = np.linspace(200,700,6)
cp = ax.contour(XG,YG,ZG,cval)
ax.clabel(cp, inline=True, fontsize=10, fmt="%d")

# Western contact
pi = np.pi
strike, dip = np.radians([20, 22]) 
point1 = np.array([692, 1212, 600])
DG = outcrop_trace(strike,dip,point1,XG,YG,ZG)
cval = 0 # Contour only CG zero value
cp = ax.contour(XG,YG,DG,cval,colors="red",linewidths=3)
ax.plot(point1[0],point1[1],"bo",markersize=10)

# Eastern contact
strike, dip = np.radians([160, 22]) 
point2 = np.array([3203, 1031, 200])
DG = outcrop_trace(strike,dip,point2,XG,YG,ZG)
cp = ax.contour(XG,YG,DG,cval,colors="red",linewidths=3)
ax.plot(point2[0],point2[1],"bo",markersize=10)

# Make axes equal
ax.axis("equal")
plt.show()


# In[ ]:




