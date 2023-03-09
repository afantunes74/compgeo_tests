#!/usr/bin/env python
# coding: utf-8

# (ch03_nb02)=
# # Maximum strike error vs. dip graph

# In[1]:


# Import libraries
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


# Make a figure
fig, ax = plt.subplots()

# Define dip
dip = np.radians(np.arange(0.1,40,0.1))
# Define operator error
eo=np.radians(np.arange(1,6,1))
for i in range(eo.size):
    # maximum strike error
    es = np.tan(eo[i])/np.tan(dip)
    # Plot each curve
    ax.plot(dip*(180/np.pi),es*(180/np.pi),"k")

# One example
# dip = 5 deg
dip_ex = 5 
# maximum operator error = 2 deg
eo_ex = 2 
# maximum strike error
pi = np.pi
es_ex = np.tan(eo_ex*pi/180)/np.tan(dip_ex*pi/180)*180/pi
# Plot example
ax.plot([dip_ex,dip_ex,0], [0,es_ex,es_ex], "k--")
ax.plot(dip_ex,es_ex,"ko")
# Add labels to eo curves
ax.text(8,7.5,"1")
ax.text(9,13,"2")
ax.text(10,17.5,"3")
ax.text(11,21,"4")
ax.text(12,24,"5")
ax.text(33,9,r"$\epsilon_o(\circ)$")
# Axes limits
ax.axis ([0, 40, 0, 90])
# Label axes
ax.set_xlabel(r"dip $\delta(\circ)$")
ax.set_ylabel(r"Maximum strike error $\epsilon_s(\circ)$")
# Add grid
ax.grid()
plt.show()


# In[ ]:




