#!/usr/bin/env python
# coding: utf-8

# (ch02_nb05)=
# # Transform locations from lat-long to World Mercator

# We have a csv file with two columns: longitude and latitude. Each coordinate pair is the center of a volcano around the world. There are 1,509 volcanoes in our dataset. The original coordinate reference system is geographic coordinates with datum WGS84. We want to make a coordinate transformation of these data points to World Mercator. It will take much too long to manually transform these coordinates as we have done in the notebooks before. Therefore, our new code will read the csv file and create a new csv file.
# 
# Check that the input (`src_dir`) and output (`dst_dir`) directories match the directory where the csv file is. In this example, the volcanoes file `volc_longlat.csv`) is in the directory `data/ch2-5`. Run the code, you will know the process is finished when the message "process completed" and the time of execution are returned:

# In[1]:


# Thanks to Rustam Zaitov for implementing 
# this new version of the code

# Import libraries
import csv, time
from os import path
from pyproj import Transformer, CRS


# In[2]:


src_file = "volc_longlat.csv" # input file
dst_file = "volc_projected.csv" # output file

src_dir = path.abspath("../../data/ch2-5") # input directory
dst_dir = path.abspath("../../data/ch2-5") # output directory

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

print("process completed in {} seconds"
      .format(end_time-start_time))


# It takes less than one second to run this code! Check the newly created csv file and notice that you now have a listing of coordinates in meters. The EPSG definition of the output coordinate reference system is listed under `dst_crs`. You can easily change this variable to another EPSG and rerun the script. If you wish to run the script on another file, change the`src_file` and `dst_file`, and the `scr_dir` and `dst_dir` if the file is in another directory.

# In[ ]:




