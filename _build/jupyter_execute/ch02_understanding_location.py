#!/usr/bin/env python
# coding: utf-8

# (ch02)=
# # Understanding location

# (ch02-1)=
# ## Location
# 
# Understanding where we are located or where an object of interest is located is very important in geosciences. Geographic information science is the study of geographic information; it includes theory and concepts and provides methods to combine and analyze spatial data (Watson, 2017). Geographic information systems is the software and technology that supports the application of geographic information theory (Watson, 2017). GIS is an acronym used for either geographic information systems or geographic information science.
# 
# Geosciences are Earth-based and so location-based. The geographic aspect is highly applicable for the geosciences. The basic usage of GIS is cartographic – making maps. In this chapter, you will become familiar with some of the basic concepts in defining locations on the Earth. These concepts are also used for defining locations on other planets as well.

# (ch02-2)=
# ## Geodesy basics
# Geodesy is “the branch of mathematics dealing with the shape and area of the Earth or large portions of it” (Lexico, 2019); however, the discipline is expanding to include other planets. We need geodesy to model the Earth and make calculations because the Earth is not a perfect sphere nor flat.

# (ch02-2-1)=
# ### Scales
# 
# In cartography, we use scales to describe how the real world length has been reduced to fit on a page or screen. Usually, map scales are described as ratios, such as 1:50,000. When we describe the scale, we say it is either small or large. The description of small or large refers to the scale, not the size of the area. Therefore, if we describe a scale as small, it means the fraction described by the scale ratio is small; this in turn means the map covers a large area. The inverse is true for large scales; the fraction described by the scale is large and in turn covers a small area (Kraak and Ormeling, 2003, Lisle et al., 2011). There is not an official categorical differentiation between small and large scales. Generally speaking, small scale maps cover regions, countries, and continents, while large scale maps cover neighborhoods, towns, or counties. The following example illustrates this.

# (ch02-2-2)=
# ### Example 1: Understanding Map Scales
# 
# A map has a scale of 1:1,000,000. Would you refer to this as a small or large scale?
# 
# 1.  First, rewrite this as a fraction: 1/1,000,000
# 2.  Is this a small fraction or a large fraction? This is a relatively small fraction, so it is a small scale.
# 
# {numref}`Figure %s <ch02_fig01>` is an example of a map with a scale of 1:1,000,000.

# ```{figure} /figures/ch02_fig01.png
# :width: 450px
# :name: ch02_fig01
# 
# A historical map of northern Scotland at the scale of 1:1,000,000 (Geographic Section - General Staff, 1941).
# ```

# A map has a scale of 1:10,000. Would you refer to this as a small or large scale?
# 
# 1.  First, rewrite this as a fraction: 1/10,000
# 2.  Is this a small fraction or a large fraction? This is a relatively large fraction, so it is a large scale.
# 
# {numref}`Figure %s <ch02_fig02>` is an example of a map with a scale of 1:10,000.

# ```{figure} /figures/ch02_fig02.png
# :width: 500px
# :name: ch02_fig02
# 
# A historical map of Rome, Italy at the scale of 1:10,000 (C.I.U. and War Office, 1944).
# ```

# (ch02-2-3)=
# ### Authalic Sphere
# 
# The authalic sphere is a sphere used for the basic surface for mapping ({numref}`Figure %s <ch02_fig03>`); its surface area is the same as the ellipsoid (Robinson et al., 1995). Based on the WGS-84 ellipsoid, the Earth has an equatorial radius of 6371 km and a circumference of 40,030.2 km. The radius is an often-used constant in geodesy (Robinson et al., 1995). The authalic sphere is used in small scale mapping (small scale covers a large area) because the difference between the authalic sphere and the ellipsoid is minimum over large areas (Robinson et al., 1995).

# ```{figure} /figures/ch02_fig03.png
# :width: 450px
# :name: ch02_fig03
# 
# Highly stylized comparison of sphere, ellipse, and geoid.
# ```

# (ch02-2-4)=
# ### Ellipsoidal Earth
# 
# Due to gravity, the Earth flattens at the poles. In a cross-section, the Earth looks like an oblate ellipsoid ({numref}`Figure %s <ch02_fig03>`) (Robinson et al., 1995). Oblateness refers to flatness. When mapping over small areas (large scale mapping), the oblateness of the ellipsoid must be taken into account. The GPS network uses the WGS-84 ellipsoid (Robinson et al., 1995).
# 
# There are varying ellipsoidal measurements on different continents and times. These continental differences are due to gravity. Temporal differences are due to technological accuracy. The WGS-84 ellipsoid is based on satellite observations and is accepted as being highly accurate (Robinson et al., 1995).

# (ch02-2-5)=
# ### Geoid
# 
# Geoid means Earth-like and is in 3D. It is based on an equipotential gravity surface. The geoid follows the mean sea level in oceans and hypothetical sea-level canals on the continents (Robinson et al., 1995). Due to geology (rock density) and topography, the geoid deviates from the ellipsoid ({numref}`Figure %s <ch02_fig01>`). The geoid is a “reference surface for ground surveyed horizontal and vertical positions” (Robinson et al., 1995).

# (ch02-3)=
# ## Projections
# 
# A projection is a mathematical equation to transfer a region, of whatever size, of the round Earth onto a flat surface. Projections are used because distance and surface area calculations are more difficult on a sphere. A flat map can show greater detail than a sphere and is more transportable. Imagine how large a globe you would need to sufficiently show the streets in your neighborhood! We need projections to transform our 3D ellipsoidal Earth onto a flat map. Projections may be based on the authalic sphere, ellipsoid, or geoid.
# 
# Before proceeding, take a moment to look over an informational [pictographic](https://pubs.usgs.gov/gip/70047422/report.pdf) by the U.S. Geological Survey describing different types of projection.

# (ch02-3-1)=
# ### Distortions
# 
# All projections have distortions that vary by projection type (i.e. transverse Mercator vs. Miller cylindrical – see pictograph mentioned in previous section). Selecting a projection depends on discipline, size of area, orientation of area, regional standards, map purpose, and map scale. There are many resources for determining which projection you should use. Large-scale mapping uses conformal projection because angles measured on the ground are the same as those in the map (Iliffe and Lott, 2008). Four types of distortion are: area, shape, direction, and distance. The Tissot’s Indicatrix is a graphic device to show the distortion at a point (Robinson et al., 1995). We will investigate this phenomenon using the Python library [Cartopy](https://scitools.org.uk/cartopy/docs/latest/#). To install this library, follow the steps below:

# ```{admonition} Installing Cartopy
# We will use the Cartopy library to visualize projection distortions using the Tissot’s Indicatrix. We will make a special Cartopy Environment in Anaconda. This is because Cartopy dependencies are lower than some other libraries we’ll be using. This happens from time to time when we use open-source libraries.
# 
# 1.  Open Anaconda Navigator.
# 2.  In the left panel, click on `Environments`.
# 3.  In the middle panel, click `Create` to create a new environment
# 4.  Name this environment: `ch2cartopy`
# 5.  Select Python version 3.8 or later
# 6.  Click `Create`. This will take a few minutes.
# 7.  We need to install Cartopy. In the right panel in the search field, type: `cartopy`
# 8.  You will get the message `0 packages available matching cartopy` since Cartopy is not installed.
# 9.  In the right panel, change `Installed` to `Not installed`
# 10. Cartopy now shows up. Click the check box next to `cartopy`.
#     
#     :::{note}
#     If it does not, press `Update index` in the right panel.
#     :::
# 
# 11. In the screen lower right corner, click `Apply`
# 12. Please wait while Cartopy is collected and then click `Apply`
# 13. The Cartopy library and any dependencies will be installed or updated. This may take a few minutes.
#     
#     :::{note}
#     If Cartopy is lower than v. 0.21, you will need to downgrade the Matplotlib from 3.6.2 to 3.5.2. Open a terminal or command window. Type `conda activate ch2cartopy` followed by `Enter`, and then `conda install matplotlib==3.5.2` followed by `Enter`.
#     :::
# 
# 14. Click on `Home`
# 15. Click on `Install` under Jupyter Notebook.
# 16. Click on `Launch` under Jupyter Notebook.
# ```

# (ch02-3-2)=
# ### Example 2: Tissot’s Indicatrix
# 
# This example will introduce you to understanding distortions in projections. As you change the projection name, a different mathematical equation will be used to portray the round Earth in a flat presentation. Pay particular attention to the size, shape, and spacing of the ellipses describing the distortion. The Tissot’s Indicatrix quickly and easily visualizes the changes in area and spatial relationships between different projections.
# 
# The notebook [ch2-1](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch2-1.ipynb) contains this example. If you just launched Jupyter Notebook as indicated above, open the notebook. If you closed Anaconda, follow these steps:
# 
# 1.  Open Anaconda Navigator
# 2.  Click on `Environments`
# 3.  Choose `ch2cartopy`
# 4.  Click `Home`
# 5.  Launch Jupyter Notebook
# 6.  Open the notebook `ch2-1`
# 
# This example starts with the Plate Carree projection ([Cartopy, 2018b](https://scitools.org.uk/cartopy/docs/latest/gallery/tissot.html)). Run the code below:

# In[1]:


# Import cartopy and matplotlib
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

# Hide warnings
import warnings 
warnings.simplefilter("ignore")


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


# The Tissot’s Indicatrix is symbolized by the orange ellipses. Closer to the poles, the ellipses become more oblate; while closer to the Equator, they are more circular. The Plate Carrée projection is a specific form of the Equidistant Cylindrical projection. Plate Carrée has the latitude of origin at the Equator.
# 
# Change the code to use the Mollweide projection:

# In[2]:


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


# Notice how the sizes of the ellipses are very similar throughout the map. All of the ellipses are more rounded and circular regardless of their position, as compared to the Plate Carrée projection. Mollweide is often used for world maps.
# 
# Cartopy has a list of projections that are included in the library ([Cartopy, 2018a](https://scitools.org.uk/cartopy/docs/latest/crs/projections.html)). Change the code above to project the map in Azimuthal Equidistant.

# (ch02-4)=
# ## Reference systems and datums
# 
# A coordinate reference system is a coordinate system that has been referenced to a datum. A datum is the location used for a reference point from which spatial measurements are made. There are geographic and Cartesian coordinate systems. Coordinates are for specific locations on the Earth. They can be expressed as geographic using latitude and longitude. Latitude are parallels that are evenly spaced and longitude are meridians that converge at the poles. These are measured in degrees. Cartesian coordinates are expressed in x and y and may have units that are meters, feet, or kilometers, for example. Coordinates only have meaning when the coordinate system and datum are known (Iliffe and Lott, 2008, Robinson et al., 1995).

# (ch02-4-1)=
# ### Example 3: Defining the coordinate reference system
# 
# In any GIS program, including spatial libraries and code, the user must ensure that the coordinate reference system is defined for the spatial data. The GIS software will make assumptions, sometimes erroneous, if the coordinate reference system is not properly defined. In this example, we will see a demonstration of these assumptions and how to prevent them. This example is from the SciTools tutorials to understand Cartopy ([Cartopy, 2018c](https://scitools.org.uk/cartopy/docs/latest/tutorials/understanding_transform.html)).
# 
# In Cartopy, there are two keywords that you must understand in order to properly display your data. The “projection” argument is used for display of your data. This only affects the map or plot. It does not define the coordinate reference system of the data itself. The “transform” argument, on the other hand, defines the coordinate reference system. The best practice is to define both of these. We will investigate the error that occurs when the best practice is not followed and compare this to when the best practice is followed. The notebook [ch2-2](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch2-2.ipynb) contains this example.
# 
# First we will create some dummy data on a regular latitude/longitude grid:
# 
# :::{note}
# In Python, the backslash character `\` can be used to split a long line of code.
# :::

# In[3]:


import numpy as np


lon = np.linspace(-80, 80, 25)
lat = np.linspace(30, 70, 25)
lon2d, lat2d = np.meshgrid(lon, lat)

data = np.cos(np.deg2rad(lat2d) * 4) + \
    np.sin(np.deg2rad(lon2d) * 4)


# In order to demonstrate the error before "best practice", we will create a map using the Plate Carree projection but only specify the `projection` argument. Remember that the best practice requires both `projection` and `transform` arguments to be defined.

# In[4]:


# Import cartopy and matplotlib
import cartopy.crs as ccrs
import matplotlib.pyplot as plt

# Hide warnings
import warnings 
warnings.simplefilter("ignore")


# The projection keyword determines how the plot will look
fig = plt.figure(figsize=(6, 3))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_global()
ax.coastlines()

# didn't use transform, but looks ok...
ax.contourf(lon, lat, data)


# In this case, the data just happen to fall in the correct location. Now, we will define the data coordinate reference system (first line of code below) and add the `transform` argument to the plot (second last line of code below).

# In[5]:


# The data are defined in lat/lon coordinate system, 
# so PlateCarree() is the appropriate choice:
data_crs = ccrs.PlateCarree()

# The projection keyword determines how the plot will look
fig = plt.figure(figsize=(6, 3))
ax = plt.axes(projection=ccrs.PlateCarree())
ax.set_global()
ax.coastlines()

# use transform
ax.contourf(lon, lat, data, transform=data_crs)


# You will notice that our plot remains unchanged. The assumption, as stated previously, is that if the coordinate reference system of the data is undefined, it is the same as the map (or plot) projection. In the example above, this has been the case. Let us now investigate what happens when changing the projection of the map without defining the coordinate reference system of the data. We now define the projection to `RotatedPole` and omit the `transform` argument to see what happens:

# In[6]:


# Now we plot a rotated pole projection
projection = ccrs.RotatedPole(pole_longitude=-177.5, pole_latitude=37.5)
fig = plt.figure(figsize=(6, 3))
ax = plt.axes(projection=projection)
ax.set_global()
ax.coastlines()

# didn't use transform, uh oh!
ax.contourf(lon, lat, data)


# In this case, we see that the country boundaries have rotated and changed shape, however, the data did not move with the the country boundaries. We made a wrong assumption in the data definition. Therefore we need to define the `transform` argument:

# In[7]:


# A rotated pole projection again...
projection = ccrs.RotatedPole(pole_longitude=-177.5, pole_latitude=37.5)
fig = plt.figure(figsize=(6, 3))
ax = plt.axes(projection=projection)
ax.set_global()
ax.coastlines()

# ...but now using the transform argument
ax.contourf(lon, lat, data, transform=data_crs)


# Now the data are correctly projected. Get in the habit of always defining the coordinate reference system and the map plot projection. It will save headaches and misunderstandings of the data.
# 
# Here is another script using an entirely different projection and additional plotting parameters. Can you figure out what the script will produce before running the cell? Notice the presence of the Matplotlib function `subplot`. What do you think this function does?

# In[8]:


# We can choose any projection we like...
projection = ccrs.InterruptedGoodeHomolosine()
fig = plt.figure(figsize=(6, 7))
ax1 = plt.subplot(211, projection=projection)
ax1.set_global()
ax1.coastlines()
ax2 = plt.subplot(212, projection=ccrs.NorthPolarStereo())
ax2.set_extent([-180, 180, 20, 90], crs=ccrs.PlateCarree())
ax2.coastlines()

# ...as long as we provide the correct transform, 
# the plot will be correct
ax1.contourf(lon, lat, data, transform=data_crs)
ax2.contourf(lon, lat, data, transform=data_crs)


# Now that you have worked through the exercise, you should have an understanding of why it is important to fully define the coordinate reference system and the projection when creating a map.

# (ch02-5)=
# ## Conversion vs. transformation
# 
# When working with data collected from several sources or in different coordinate reference systems, the data must be redefined to have the same coordinate reference system and datum. Consistent coordinate reference systems for data in a map is important because there may be spatial differences between the coordinate reference systems creating locational errors. A coordinate conversion is when the coordinate reference systems have the same datum. A coordinate transformation is when the coordinate reference systems have different datums.
# 
# When making a map, all of the data in the map should have the same coordinate reference system definition. Software include many definitions and transformations but refer to Snyder (1987) and International Association of Oil and Gas Producers (2018) for projection formulae (Iliffe and Lott, 2008). In Example 2, the Tissot’s Indicatrix ellipses did not change coordinate reference systems. The method used to plot the ellipses changed. In Example 3, we also only changed the projection of the map, not the coordinate reference system. In the following examples, we will make coordinate conversions and transformations. For this purpose, we will use the [pyproj](https://pypi.org/project/pyproj/) Python library and examples from the pyproj documentation ([Whitaker, 2019](https://pyproj4.github.io/pyproj/stable/api/transformer.html?highlight=transformer)).

# ```{admonition} Installing pyproj
# 
# A version of pyproj was installed with Cartopy, but you need the latest version which has different dependencies than Cartopy. You will create a new environment.
# 
# 1.  Open Anaconda Navigator.
# 2.  In the left panel, click on `Environments`
# 3.  In the middle panel, click `Create` to create a new environment
# 4.  Name this environment: `ch2pyproj`
# 5.  Select Python version 3.8 or later
# 6.  Click `Create`. This will take a few minutes
# 7.  We need to install pyproj. In the right panel in the search field, type: `pyproj`
# 8.  You will get the message `0 packages available matching pyproj` since pyproj is not installed.
# 9.  In the right panel, change `Installed` to `Not installed`
# 10. pyproj now shows up. Click the check box next to `pyproj`.
# 11. In the screen lower right corner, click `Apply`
# 12. Please wait while pyproj is collected and then click `Apply`
# 13. The pyproj library and any dependencies will be installed or updated. This may take a few minutes.
# 14. Click on `Home`
# 15. Click on `Install` under Jupyter Notebook.
# 16. Click on `Launch` under Jupyter Notebook.
# ```

# (ch02-5-1)=
# ### Example 4: Coordinate conversion
# 
# The notebook [ch2-3](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch2-3.ipynb) contains this example. If you just launched Jupyter Notebook as indicated above, open the notebook. If you closed Anaconda, follow these steps:
# 
# 1.  Open Anaconda Navigator
# 2.  Click on `Environments`
# 3.  Choose `ch2pyproj`
# 4.  Click `Home`
# 5.  Launch Jupyter Notebook
# 6.  Open the notebook `ch2-3`
# 
# We have a coordinate pair defined in decimal degrees of latitude and longitude. The longitude is -120.108and latitude is 34.36116666. We want to make a coordinate conversion from latitude and longitude to Universal Transverse Mercator, where the point is defined by east and north coordinates in meters. To learn more about Universal Transverse Mercator (UTM), refer to (Snyder, 1987). In the code, we use the pyproj `Proj` function. We can only use `Proj` when making a coordinate conversion (i.e. the same datum):

# In[9]:


# Import pyproj
from pyproj import Proj


# Construct the projection matrix
p = Proj(proj="utm",zone=10,ellps="WGS84", 
         preserve_units=False)

# Apply the projection to the lat-long point
x,y = p(-120.108, 34.36116666)

print(f"x={x:9.3f}, y={y:11.3f}")


# This is the same location but only expressed in east and north coordinates in meters using the UTM coordinate reference system. The datum used is WGS84. We can convert the UTM coordinates back to latitute and longitude by adding two lines of code:

# In[10]:


# Apply the inverse of the projection matrix
# to the point in UTM
lon,lat = p(x,y,inverse=True)
print(f"lon={lon:8.8f}, lat={lat:5.8f}")


# We can confirm that the inverse conversion arrives at the original pair. Let’s now try converting several points of different latitude and longitude using a collection of objects in Python, or tuples. Add the following code:

# In[11]:


# three points in lat-long
lons = (-119.72,-118.40,-122.38)
lats = (36.77, 33.93, 37.62 )
# Apply the projection to the points
x1,y1 = p(lons, lats)
print(x1,y1)


# Now, let’s do a more advanced exercise: In the cartographic community, an easy way to communicate the coordinate reference system is to use the EPSG Geodetic Parameter Data set. Every coordinate reference system is given a code. This ensures that if someone uses UTM zone 10 North with datum WGS-84 and tells you UTM zone 10, that you do not accidentally use UTM zone 10 North with datum GRS80, for example.
# 
# Earlier in this exercise, we defined the UTM zone in the `Proj` function. Here, we will refer to the EPSG code. First, we will take a coordinate pair in longitude and latitude with datum WGS84 and convert it to EPSG:32667. Before proceeding, conduct a quick internet search on what EPSG:32667 means. This is important to understand what we will do next. The first part of the code is:

# In[12]:


# silence warnings
import warnings 
warnings.simplefilter("ignore")


# initial coordinate conversion
p = Proj(init="EPSG:32667", preserve_units=True, 
         always_xy=True)
# Apply the conversion to the lat-long point
x,y = p(-114.057222, 51.045)
print(f"x={x:9.3f}, y={y:11.3f} (feet)")


# Let’s dissect this as the pyproj code looks quite a bit different. The first part of the function `Proj` calls EPSG:32667. If you looked up EPSG:32667 online, you found that it is for UTM zone 17 North, but the units are in feet. The default mode for `Proj` is `preserve_units=False`, which forces any unit to meters. However, we want to see the units in US Survey Feet as the projection defines; therefore, we change the argument to `True`.
# 
# Now, suppose we want to see the output in meters. How will you amend the code? Here is what you should add:

# In[13]:


# Print the coordinate pair in meters
p1 = Proj(init="EPSG:32667", preserve_units=False)
x1,y1 = p1(-114.057222, 51.045)
print(f"x={x1:9.3f}, y={y1:11.3f} (meters)")


# As discussed, you should change `preserve_units=False` and change the unit to be printed from `feet` to `meters`. Congratulations! You now have a good understanding of coordinate conversions.

# (ch02-5-2)=
# ### Example 5: Coordinate transformation
# 
# We learned earlier that we have a coordinate conversion where a coordinate pair is converted between coordinate reference systems with the same datum. In many instances, the coordinate reference system will also undergo a datum shift – this is a coordinate transformation.
# 
# This example is included in the notebook [ch2-4](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch2-4.ipynb). We will use the pyproj `CRS` and `transform` functions. The `CRS` function defines the coordinate reference system while the `transform` function specifies which coordinate reference system is the original and which is the output. `CRS` has the same ability to refer directly to an EPSG code.
# 
# The input coordinates are in EPSG:4326, which is a commonly used code. It is the geographic coordinate system with datum WGS84. The output coordinates are EPSG:31984, which is for UTM zone 24 S with datum SIRGAS2000.

# In[14]:


# Import transform and CRS functions
from pyproj import transform
from pyproj import CRS

# Function transform is deprecated
# Silence warnings
import warnings 
warnings.simplefilter("ignore")


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


# (ch02-5-3)=
# ### Example 6: Transforming several points at once
# 
# We have focused our examples on one coordinate pair at a time. The reality is that you will more often have several coordinates to transform at one time. The notebook [ch2-5](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch2-5.ipynb) explains how to do this.
# 
# We have a csv file with two columns: longitude and latitude. Each coordinate pair is the center of a volcano around the world. There are 1,509 volcanoes in our dataset. The original coordinate reference system is geographic coordinates with datum WGS84. We want to make a coordinate transformation of these data points to World Mercator. It will take much too long to manually transform these coordinates as we have done in the notebooks before. Therefore, our new code will read the csv file and create a new csv file.
# 
# :::{note}
# Notice that the code can be further simplified by the use of the Numpy library or the Pandas library.
# :::
# 
# Check that the input (`src_dir`) and output (`dst_dir`) directories match the directory where the csv file is. In this example, the volcanoes file ([volc_longlat.csv](https://github.com/nfcd/compGeo/blob/master/source/data/ch2-5/volc_longlat.csv)) is in the directory `data/ch2-5`. Run the code, you will know the process is finished when the message `process completed` and the time of execution are returned:

# In[15]:


# Thanks to Rustam Zaitov for implementing 
# this new version of the code

# Import libraries
import csv, time
from os import path
from pyproj import Transformer, CRS


src_file = "volc_longlat.csv" # input file
dst_file = "volc_projected.csv" # output file

src_dir = path.abspath("data/ch2-5/") # input directory
dst_dir = path.abspath("data/ch2-5/") # output directory

src_path = path.join(src_dir, src_file)
dst_path = path.join(dst_dir, dst_file)

src_crs = CRS("EPSG:4326") #WGS84
dst_crs = CRS("EPSG:3395") #World Mercator

# create coordinate transformer
# always_xy=True makes projector.transform() accept 
# lon, lat (GIS order) instead of lat, lon
projector = Transformer.from_crs(src_crs, dst_crs, 
                                 always_xy=True)

# source csv file has lon, lat columns
src_header = ["LONGITUDE", "LATITUDE"]

# destinatin csv file will have x, y columns
dst_header = ["x", "y"]

# start benchmark timer
start_time = time.time()

# open destination file in write mode
with open(dst_path, "w") as w:
    # open source file in read mode
    with open(src_path, "r") as r:
        reader = csv.reader(r, dialect="excel")
        # read and skip first header row 
        input_headers = next(reader)         

        writer = csv.writer(w, delimiter=",", quotechar='"',
                            quoting=csv.QUOTE_MINIMAL)
        # Write the output header
        writer.writerow(dst_header)   
        for row in reader:
            try:
                # convert string values inside row 
                # into float values
                lon, lat = [float(val) for val in row]
                x, y = projector.transform(lon, lat)
                writer.writerow([ x, y ])
            except Exception as e:
                # If coordinates are out of bounds, 
                # skip row and print the error
                print (e)

# stop benchmarking
end_time = time.time()

print("process completed in {} seconds".format(end_time-start_time))


# It takes less than one second to run this code! Check the newly created csv file and notice that you now have a listing of coordinates in meters. The EPSG definition of the output coordinate reference system is listed under `dst_crs`. You can easily change this variable to another EPSG and rerun the script. If you wish to run the script on another file, change the `src_file` and `dst_file`, and the `src_dir` and `dst_dir` if the file is in another directory.

# (ch02-6)=
# ## Exercises
# 
# 1. So far, you have worked with data sets that have been provided to you. It is an important skill to be able to find data sets online and to prepare them for use in your work.
# 
#     The United States Geological Survey (USGS) Earthquake Hazards Program monitors, records, and maintains a global database of earthquake activity. The public can query the archive and download earthquake epicenter data. Go to the [Search Earthquake Catalog](https://earthquake.usgs.gov/earthquakes/search/) page of the USGS. Using the Basic Options, download all of the magnitude 4.5+ earthquakes in the last year in the world. In the Output Options, choose a CSV format.
# 
#     1. Modify the notebook [ch2-5](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch2-5.ipynb) to transform the earthquake epicenters from latitude-longitude to World Mercator (or another map projection of your choice). Use this [pictographic](https://pubs.usgs.gov/gip/70047422/report.pdf) for further information. The datum of the earthquake epicenters is likely WGS84.
# 
#     2. Plot the earthquake epicenters. Make sure to include the outline of the continents.
#     
#     :::{hint}
#     Look at the notebook [ch2-2](https://github.com/nfcd/compGeo/blob/master/source/notebooks/ch2-2.ipynb) for a starting point, and check the [Cartopy](https://scitools.org.uk/cartopy/docs/latest/#) website for examples of how to plot localities on a map.
#     :::
# 
#     3. Now modify your plot to color the earthquakes by depth (red are shallow and blue are deep earthquakes), and the size of the points by the earthquake magnitude.
# 
#     4. Add the volcanoes from Example 6 ([volc_longlat.csv](https://github.com/nfcd/compGeo/blob/master/source/data/ch2-5/volc_longlat.csv)) to the map (use triangles to indicate volcanoes). Do you see any correlation between the earthquake epicenters and the volcanoes?

# (ch02-7)=
# ## References
# 
# C.I.U. AND WAR OFFICE. 1944. Town Plan of Roma (Rome) (North Sheet),
# 1:10,000. Washington, D. C.: War Office.
# 
# Cartopy. 2018a. Projections
# \[[Online](https://scitools.org.uk/cartopy/docs/latest/crs/projections.html)\].
# UK: SciTools. \[Accessed 19 November, 2019\]
# 
# Cartopy. 2018b. Tissot’s Indicatrix
# \[[Online](https://scitools.org.uk/cartopy/docs/latest/gallery/tissot.html)\].
# UK: SciTools. \[Accessed 19 November, 2019\].
# 
# Cartopy. 2018c. Understanding the Transform and Projection Keywords
# \[[Online](https://scitools.org.uk/cartopy/docs/latest/tutorials/understanding_transform.html)\].
# UK: SciTools. \[Accessed 19 November, 2019\].
# 
# Geographic Section - General Staff. 1941. Aberdeen, 1:1,000,000. Great
# Britain: War Office.
# 
# Iliffe, J. and Lott, R. 2008. Datums and Map Projections: For Remote
# Sensing, GIS, and Surveying, Dunbeath, Scotland, Whittles.
# 
# International Association of Oil and Gas Producers. 2018. Geomatics
# Guidance Note 7, Part 2 Coordinate conversions and Transformations
# including Formulas.
# 
# Kraak, M.J. and Ormeling, F.J. 2003. Cartography: visualization of
# geospatial data. Addison Wesley.
# 
# Lexico. 2019. Geodesy
# \[[Online](https://www.lexico.com/en/definition/geodesy)\]. Oxford.
# \[Accessed August, 2019\].
# 
# Lisle, R. J., Brabham, P. and Barnes, J. W. 2011. Basic Geological
# Mapping, Chicester, UK, Wiley-Blackwell.
# 
# Robinson, A. H., Morrison, J. L., Muehrcke, P. C., Kimerling, A. J. and
# Guptill, S. C. 1995. Elements of Cartography, New York, Wiley.
# 
# Snyder, J. P. 1987. Map Projections: a working manual. Geological Survey
# Professional Paper. Washington, D. C., U.S.A.: United States Government
# Printing Office.
# 
# Watson, L. 2017. Spatial-based assessment at continental to global
# scale: case studies in petroleum exploration and ecosystem services.
# PhD, Utrecht University.
# 
# Whitaker, J. 2019. pyproj Transformer Documentation
# \[[Online](https://pyproj4.github.io/pyproj/stable/api/transformer.html?highlight=transformer)\].
# \[Accessed 7 January, 2020\].

# In[ ]:




