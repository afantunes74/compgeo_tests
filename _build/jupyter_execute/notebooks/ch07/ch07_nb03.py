#!/usr/bin/env python
# coding: utf-8

# (ch07_nb03)=
# # The Mohr circle for stress in 3D

# In[1]:


# Import libraries
import numpy as np
pi = np.pi

# Import mohr_circle_stress
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.mohr_circle_stress import mohr_circle_stress


# In[2]:


# Stress tensor in principal stress coordinate system
stress = np.array([[50, 0, 0],[ 0, 30, 0],[ 0, 0, 10]])

# Trend and plunge of sigma1, and trend of sigma3
tx1, px1, tx3 = np.radians([90, 0, 90])

# Planes
planes = np.zeros((12,2))
# Strikes in degrees
planes[0:3,0] = 0
planes[3:6,0] = 180
planes[6:9,0] = 45
planes[9:12,0] = 135
# Dips in degrees
planes[0:12:3,1] = 30
planes[1:12:3,1] = 45
planes[2:12:3,1] = 60

# Convert to radians
planes = planes * pi/180
 
# Plot Mohr circle
ns,ons,fig,ax=mohr_circle_stress(stress,tx1,px1,tx3,planes)

# Print normal and shear tractions
print("Strike","Dip","\u03C3","Trend","Plunge","\u03C4",
      "Trend","Plunge",sep="\t")

# return to degrees
planes = planes*180/pi
ons = ons*180/pi
# print
for i in range(0,np.size(planes,0)):
    print(f"{planes[i,0]:.1f}",f"{planes[i,1]:.1f}",
          f"{ns[i,0]:.2f}",f"{ons[i,0]:.1f}",
          f"{ons[i,1]:.1f}",f"{ns[i,1]:.2f}",
          f"{ons[i,2]:.1f}",f"{ons[i,3]:.1f}",sep="\t")


# In[3]:


# Run this cell if you want to save the figure
fig.savefig("mohr_circle_stress.png", dpi=300)


# In[ ]:




