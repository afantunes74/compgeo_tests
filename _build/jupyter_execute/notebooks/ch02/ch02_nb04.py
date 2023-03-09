#!/usr/bin/env python
# coding: utf-8

# (ch02_nb04)=
# # Coordinate transformation

# We will use the pyproj `CRS` and `transform` functions. The `CRS` function defines the coordinate reference system while`transform` specifies which coordinate reference system is the original and which is the output.`CRS` has the same ability to refer directly to an EPSG code.
# 
# The input coordinates are in EPSG:4326, which is a commonly used code. It is the geographic coordinate system with datum WGS84. The output coordinates are EPSG:31984, which is for UTM zone 24 S with datum SIRGAS2000.

# In[1]:


# Import transform and CRS functions
from pyproj import transform
from pyproj import CRS

# Function transform is deprecated
# Silence warnings
import warnings 
warnings.simplefilter("ignore")


# In[2]:


# input coordinates
c1 = CRS("EPSG:4326")
# coordinate pair
y1=-10.754283
x1=-39.866132
# output coordinates
c2 = CRS("EPSG:31984")
# Coordinate transformation
x2, y2 = transform(c1, c2, x1, y1)
print(f"x={x2:9.3f}, y={y2:11.3f}")


# In[ ]:




