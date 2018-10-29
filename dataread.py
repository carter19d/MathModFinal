
#From https://www.earthdatascience.org/courses/earth-analytics-python/lidar-raster-data/open-lidar-raster-python/
import rasterio as rio
from rasterio.merge import merge
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

# Open raster data
#lidar_dem = rio.open('n42w112_10m.tif')
lidar_dem1 = rio.open('USGS_NED_1_n41w112_IMG.img')
lidar_dem2 = rio.open('USGS_NED_1_n42w112_IMG.img')
a = lidar_dem1.read(1)
b = lidar_dem2.read(1)
print(a.shape, b.shape)

import math

pi = math.pi
cos = lambda x: math.cos(x/180*pi)
acos = lambda x: 180/pi*math.acos(x)
sin = lambda x: math.sin(x/180*pi)
asin = lambda x: 180/pi*math.asin(x)
tan = lambda x: math.tan(x/180*pi)
atan = lambda x: 180/pi*math.atan(x)

def slope_at_point(mat, point, length=10):
    try:
        el = mat[point]
    except:
        print(point)
    maxslope = 0
    maxdir = 0
    for direction in range(0, 360, 5):
        newpoint = (round(point[0] + length*cos(direction)), round(point[1] + length*sin(direction)))
        newel = mat[newpoint]
        #print(newel-el)
        slope = (newel-el)/(length)
        if slope > maxslope:
            maxslope = slope
            maxdir = direction
    #print(point, maxslope, maxdir)
    return (maxslope, maxdir)

def angle_for_slope(mat, point, desired, curdirection, scale=23.4, length=5):
    '''el = mat[point]
    bestdev = 90+desired
    bestpoint = None
    bestdir = None
    #for direction in range(curdirection - 45, curdirection + 45):
    for direction in range(0, 360):
        newpoint = (round(point[0] + length*cos(direction)), round(point[1] + length*sin(direction)))
        newel = mat[newpoint]
        slope = (el-newel)/length
        if abs(atan(slope)-desired) < bestdev and atan(slope) > -5:
            bestdev = abs(atan(slope)-desired)
            bestpoint = newpoint
            bestdir = direction'''
    r = slope_at_point(mat, point, length)
    if tan(desired) <= r[0]:
        direction = r[1] + acos(tan(desired)/r[0])
        if mat[round(point[0] + length*cos(direction)), round(point[1] + length*sin(direction))] > mat[point]:
            direction = - r[1] - acos(tan(desired)/r[0])
    else:
        direction = r[1]
    return ((round(point[0] + length*cos(direction)), round(point[1] + length*sin(direction))), direction)
    #return (bestdev, bestpoint, bestdir)

def eulers(mat, point, desired, n=100, **kwargs):
    points = [point]
    curdirection = 0
    while True:
        r = angle_for_slope(mat, point, desired, curdirection, **kwargs)
        if r is None:
            return points
        if curdirection == r[1]:
            break
        curdirection = r[1]
    for _ in range(n):
        r = angle_for_slope(mat, point, desired, curdirection, **kwargs)
        #print(r[0], slope_at_point(mat, point, length=10)[0], slope_at_point(mat, point, length=20)[0])
        if r is None:
            return points
        point = r[0]
        points.append(point)
    return points

mat = b
fig, ax = plt.subplots()
ax.matshow(b)
for y in [250, 300]:
    for x in range(1240, 1280, 5):
        point = (x, y)
        print(slope_at_point(mat, point, length=10)[0], slope_at_point(mat, point, length=20)[0])
        z = eulers(b, point, 50, length=2)
        Z = list(zip(*z))
        ax.plot(Z[0], Z[1], 'r-')
plt.show()
