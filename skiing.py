
#From https://www.earthdatascience.org/courses/earth-analytics-python/lidar-raster-data/open-lidar-raster-python/
import random
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
plt.set_cmap('pink')

# Open raster data
#lidar_dem = rio.open('n42w112_10m.tif')
lidar_dem1 = rio.open('USGS_NED_1_n41w112_IMG.img')
lidar_dem2 = rio.open('USGS_NED_1_n42w112_IMG.img')
L = [(lidar_dem1, lidar_dem1.read(1)), (lidar_dem2, lidar_dem2.read(1))]
#a = lidar_dem1.read(1)
#b = lidar_dem2.read(1)
#print(a.shape, b.shape)

import math

pi = math.pi
cos = lambda x: math.cos(x/180*pi)
acos = lambda x: 180/pi*math.acos(x)
sin = lambda x: math.sin(x/180*pi)
asin = lambda x: 180/pi*math.asin(x)
tan = lambda x: math.tan(x/180*pi)
atan = lambda x: 180/pi*math.atan(x)
atan2 = lambda y,x: 180/pi*math.atan2(y,x)
sqrt = math.sqrt

R = 6371000

def inbounds(coord, bounds):
    return bounds[0]<coord[1]<bounds[2] and bounds[1]<coord[0]<bounds[3]

def elevation_at_coord(coord, data):
    for i in data:
        if inbounds(coord, i[0].bounds):
            return i[1][i[0].index(*coord[::-1])]

def dist(coord1, coord2):
    lat1, lon1 = (i/360 for i in coord1)
    lat2, lon2 = (i/360 for i in coord2)
    dlon = lon2 - lon1
    dlat = lat2 - lat1
    a = (sin(dlat/2))**2 + cos(lat1) * cos(lat2) * (sin(dlon/2))**2
    c = 2 * atan2(sqrt(a), sqrt(1-a))
    return R * c

def step_in_dir(coord, dist, direction):
    lat1, lon1 = (i/360 for i in coord)
    lat2 = asin(sin(lat1)*cos(dist/R) + cos(lat1)*sin(dist/R)*cos(direction))
    lon2 = lon1 + atan2(sin(direction)*sin(dist/R)*cos(lat1), cos(dist/R)-sin(lat1)*sin(lat2))
    return (360*lat2, 360*lon2)

def slope_points(coord1, coord2):
    return (elevation_at_coord(coord1, L)-elevation_at_coord(coord2, L)) / dist(coord1, coord2)

def maxslope(coord):
    slope = -100000
    bestangle = None
    bestpoint = None
    for direction in range(0, 360, 10):
        newcoord = step_in_dir(coord, 20, direction)
        if not elevation_at_coord(newcoord, L) is None:
            newslope = slope_points(coord, newcoord)
            if newslope > slope:
                slope = newslope
                bestangle = direction
                bestpoint = newcoord
    return slope, bestangle, bestpoint

def target_slope(coord, desired):
    r = maxslope(coord)
    if desired <= r[0]:
        direction = r[1] + acos(desired/r[0])
        newcoord = step_in_dir(coord, 20, direction)
    else:
        newcoord = r[2]
    return newcoord

def eulers(point, desired, n=50):
    points = [point]
    curdirection = 0
    for __ in range(n):
        point = target_slope(point, desired)
        points.append(point)
    return points

def plot_points(pointsL):
    mat = np.concatenate((L[1][1], L[0][1]))[1800:5400]
    plt.matshow(mat)
    for points, color in pointsL:
        Z = [L[1][0].index(*i[::-1])[::-1] for i in points]
        Z = list(zip(*[(i,j-1800) for i,j in Z]))
        plt.plot(*Z, color=color, linestyle='-')
    plt.show()

RANGES = {'b':[.06, .25], 'i':[.25, .40], 'e':[.40, 1]}
def randdesired(difficulty):
    desiredrange = RANGES[difficulty[0]]
    return random.random() * (desiredrange[1] - desiredrange[0]) + desiredrange[0]

COLORS = {'b':'g', 'i':'b', 'e':'k'}
def points_and_color(coord, difficulty):
    return (eulers(coord, randdesired(difficulty)), COLORS[difficulty[0]])

def full_ski_lift(liftcoord, n=20, N=(5, 10, 10), difficulties='bie'):
    pointsL = []
    for i in range(n):
        coord = ((liftcoord[0][0]*i + liftcoord[1][0]*(n-1-i))/(n-1),
                 (liftcoord[0][1]*i + liftcoord[1][1]*(n-1-i))/(n-1))
        for difficulty in difficulties:
            pass#difficulty = random.choice(['b', 'i', 'i', 'e', 'e'])
            pointsL.append(points_and_color(coord, difficulty))
    actualpoints = []
    for difficulty in range(len(N)):
        bests = [pointsL[difficulty]]
        bestscores = [score(pointsL[difficulty][0], difficulties[difficulty])]
        for points, color in pointsL[difficulty::len(N)]:
            if score(points, difficulties[difficulty]) < max(bestscores):
                i = bestscores.index(max(bestscores))
                bestscores.remove(bestscores[i])
                bests.remove(bests[i])
                bestscores.append(score(points, difficulties[difficulty]))
                bests.append((points,color))
            elif len(bestscores) < N[difficulty]:
                bestscores.append(score(points, difficulties[difficulty]))
                bests.append((points,color))
        actualpoints.extend(bests)
    plot_points(actualpoints)

def score(points, difficulty):
    out = 0
    for point in points:
        slope = maxslope(point)[0]
        if slope < RANGES[difficulty[0]][0]:
            out += RANGES[difficulty[0]][0]/slope
        if slope > RANGES[difficulty[0]][1]:
            out += slope/RANGES[difficulty[0]][1]
    return out

#F=Francis Peak (nearby)
F = (41 + 1/60 + 59/3600, -(111 + 50/60 + 18/3600))
Fp = (F[0], F[1]+1/60)
Fpp = (Fp[0]+5/60, Fp[1]-.5/60)
Fp = (Fp[0], Fp[1]+.5/60)
full_ski_lift([Fp, Fpp],10)
