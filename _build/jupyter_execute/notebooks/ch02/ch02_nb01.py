#!/usr/bin/env python
# coding: utf-8

# (ch02_nb01)=
# # Tissot’s Indicatrix

# This example starts with the Plate Carree projection (Cartopy, 2018b). Run the code below:

# In[1]:


# Import cartopy and matplotlib
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

# Hide warnings
import warnings 
warnings.simplefilter("ignore")


# In[2]:


# Figure
fig = plt.figure(figsize=(10, 5))
# Plate Carree projection
ax = plt.axes(projection=ccrs.PlateCarree())

# Make the map global rather than have it zoom in to
# the extents of any plotted data
ax.set_global()

# Earth image
ax.stock_img()
# Coastlines
ax.coastlines()

# Tissot"s indicatrix: Orange ellipses
ax.tissot(facecolor="orange", alpha=0.4)
plt.show()


# The Tissot’s Indicatrix is symbolized by the orange ellipses. Closer to the poles, the ellipses become more oblate; while closer to the Equator, they are more circular. The Plate Carrée projection is a specific form of the Equidistant Cylindrical projection. Plate Carrée has the latitude of origin at the Equator.
# 
# Change the code to use the Mollweide projection:

# In[3]:


# Figure
fig = plt.figure(figsize=(10, 5))
# Mollweide projection
ax = plt.axes(projection=ccrs.Mollweide())

# Make the map global rather than have it zoom in to
# the extents of any plotted data
ax.set_global()

ax.stock_img()
ax.coastlines()

ax.tissot(facecolor="orange", alpha=0.4)
plt.show()


# Notice how the sizes of the ellipses are very similar throughout the map. All of the ellipses are more rounded and circular regardless of their position, as compared to the Plate Carree projection. Mollweide is often used for world maps.
# 
# Cartopy has a list of projections that are included in the library (Cartopy, 2018a). Change the code above to project the map in Azimuthal Equidistant.

# In[4]:


# run this cell to save figure
fig.savefig("projection.png", dpi=300)


# In[ ]:




