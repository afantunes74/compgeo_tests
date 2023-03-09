#!/usr/bin/env python
# coding: utf-8

# (ch05_nb04)=
# # Down-plunge projection

# In[1]:


# Import libraries
import numpy as np
import matplotlib.pyplot as plt

# Import function down_plunge
import sys, os
sys.path.append(os.path.abspath("../../"))

from compgeo.down_plunge import down_plunge


# In[2]:


# Trend and plunge of the fold axis in radians
trend, plunge = np.radians([125, 26]) 

# Read the tops from the text files
jtc = np.loadtxt(os.path.abspath("../../data/ch5-4/jtc.txt"))
js = np.loadtxt(os.path.abspath("../../data/ch5-4/js.txt"))
kp = np.loadtxt(os.path.abspath("../../data/ch5-4/kp.txt"))

# Transform the points
jtcdp = down_plunge(jtc,trend,plunge)
jsdp = down_plunge (js,trend,plunge)
kpdp = down_plunge(kp,trend,plunge)

# Make a figure
fig, ax = plt.subplots()

# Plot the down plunge section
ax.plot(jtcdp[:,1],jtcdp[:,0],"k-")
ax.plot(jsdp[:,1],jsdp[:,0],"r-")
ax.plot(kpdp[:,1],kpdp[:,0],"b-")

# Display the section"s orientation
# Notice that the fold axis plunges SE
# Therefore the left side of the section is NE
# and the right side is SW
ax.text(-2.1e4,12e3,"NE")
ax.text(-0.6e4,12e3,"SW")

# Make axes equal
ax.axis("equal")
plt.show()


# In[ ]:




