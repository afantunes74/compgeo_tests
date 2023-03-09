#!/usr/bin/env python
# coding: utf-8

# (ch03_nb04)=
# # Maximum trend error vs. dip graph [2]

# In[1]:


# Import libraries
import numpy as np
import matplotlib.pyplot as plt


# In[2]:


# Make a figure
fig, ax = plt.subplots()

# Define dip
dip = np.radians(np.arange(0.1,90,0.1))
# Define rake
r = np.radians(np.arange(10,90,10))
# Define maximum operator error
eo = 3 * np.pi/180;
for i in range(r.size):
    k = np.tan(r[i])-np.tan(r[i]-eo)
    q = np.tan(r[i])*np.tan(r[i]-eo)
    et = (k*np.cos(dip)/(1+q*(np.cos(dip))**2))
    # Plot et versus dip'
    ax.plot(dip*180/np.pi,et*180/np.pi,'k-')
# Axes limits
ax.axis ([0, 90, 0, 11]);
# Label axes
ax.set_xlabel(r"dip $\delta(\circ)$")
ax.set_ylabel(r"Maximum trend error $\epsilon_t(\circ)$")
# Add grid
ax.grid()
plt.show()


# In[ ]:




