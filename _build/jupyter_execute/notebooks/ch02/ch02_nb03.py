#!/usr/bin/env python
# coding: utf-8

# (ch02_nb03)=
# # Coordinate conversion

# We have a coordinate pair defined in decimal degrees of latitude and longitude. The longitude is -120.108° and latitude is 34.36116666°. We want to make a coordinate conversion from latitude and longitude to Universal Transverse Mercator, where the point is defined by east and north coordinates in meters. To learn more about Universal Transverse Mercator (UTM), refer to (Snyder, 1987). In the code, we use the pyproj `Proj` function. We can only use `Proj` when making a coordinate conversion (i.e. the same datum). Run the code below.

# In[1]:


# Import pyproj
from pyproj import Proj


# In[2]:


# Construct the projection matrix
p = Proj(proj="utm",zone=10,ellps="WGS84", 
         preserve_units=False)

# Apply the projection to the lat-long point
x,y = p(-120.108, 34.36116666)

print(f"x={x:9.3f}, y={y:11.3f}")


# This is the same location but only expressed in meters (east and north) using the UTM coordinate reference system. The datum used is WGS84. We can convert the UTM coordinates back to latitute and longitude by adding two lines of code:

# In[3]:


# Apply the inverse of the projection matrix
# to the point in UTM
lon,lat = p(x,y,inverse=True)
print(f"lon={lon:8.8f}, lat={lat:5.8f}")


# We can confirm that the inverse conversion arrives at the original pair. Let's now try converting several points of different latitude and longitude using a collection of objects in Python, or tuples. Add the following code:

# In[4]:


# three points in lat-long
lons = (-119.72,-118.40,-122.38)
lats = (36.77, 33.93, 37.62 )
# Apply the projection to the points
x1,y1 = p(lons, lats)
print(x1,y1)


# Now, let's do a more advanced exercise: In the cartographic community, an easy way to communicate the coordinate reference system is to use the EPSG Geodetic Parameter Data set. Every coordinate reference system is given a code. This ensures that if someone uses UTM zone 10 North with datum WGS-84 and tells you UTM zone 10, that you do not accidentally use UTM zone 10 North with datum GRS80, for example.
# 
# Earlier in this exercise, we defined the UTM zone in the Proj function. Here, we will refer to the EPSG code. First, we will take a coordinate pair in longitude and latitude with datum WGS84 and convert it to EPSG:32667. Before proceeding, conduct a quick internet search on what EPSG:32667 means. This is important to understand what we will do next. The first part of the code is:

# In[5]:


# silence warnings
import warnings 
warnings.simplefilter("ignore")

# initial coordinate conversion
p = Proj(init="EPSG:32667", preserve_units=True, 
         always_xy=True)
# Apply the conversion to the lat-long point
x,y = p(-114.057222, 51.045)
print(f"x={x:9.3f}, y={y:11.3f} (feet)")


# Let’s dissect this as the pyproj code looks quite a bit different. The first part of the function `Proj` calls EPSG:32667. If you looked up EPSG:32667 online, you found that it is for UTM zone 17 North, but the units are in feet. The default mode for`Proj` is `preserve_units=False`, which forces any unit to meters. However, we want to see the units in US Survey Feet as the projection defines; therefore, we change the argument to `True`.
# 
# Now, suppose we want to see the output in meters. How will you amend the code? Here is what you should add:

# In[6]:


# Print the coordinate pair in meters
p1 = Proj(init="EPSG:32667", preserve_units=False)
x1,y1 = p1(-114.057222, 51.045)
print(f"x={x1:9.3f}, y={y1:11.3f} (meters)")


# As discussed, you should change to `preserve_units=False` and change the unit to be printed from `feet` to `meters`. Congratulations! You now have a good understanding of coordinate conversions.

# In[ ]:




