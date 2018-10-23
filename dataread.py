
#From https://www.earthdatascience.org/courses/earth-analytics-python/lidar-raster-data/open-lidar-raster-python/
import rasterio as rio
from rasterio.plot import show
import matplotlib.pyplot as plt
import numpy as np
import os
#plt.ion()

from shapely.geometry import Polygon, mapping
from rasterio.mask import mask
# A package created for this class that will be discussed later in this lesson
#import earthpy as et

# Set standard plot parameters for uniform plotting
plt.rcParams['figure.figsize'] = (8, 8)

# Prettier plotting with seaborn
import seaborn as sns; 
sns.set(font_scale=1.5)
sns.set_style("white")

# Open raster data
lidar_dem = rio.open('n42w112_10m.tif')

# Plot the dem using raster.io
#fig, ax = plt.subplots()
#show(lidar_dem, ax=ax)
#ax.set_axis_off()

print(list(lidar_dem.sample([(0,0)]))[0])
