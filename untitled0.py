#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 15 14:37:21 2021

@author: theone
"""
from shapely.geometry import Point,LineString,mapping,shape,MultiLineString,MultiPolygon,Polygon,MultiPoint,LinearRing
import json
import os.path
from svgpath2mpl import parse_path
import numpy as np 
import matplotlib.pyplot as plt
from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, Close
from matplotlib.path import Path
import matplotlib.patches as patches
from geopandas import GeoDataFrame
from shapely.ops import transform
from shapely.affinity import translate
import re
import textwrap
from rdp import rdp
import math

from scipy.spatial import Voronoi, voronoi_plot_2d
points = np.array([[0, 0], [0, 1], [0, 2], [1, 0], [1, 1], [1, 2],
                  [2, 0], [2, 1], [2, 2]])

border = Polygon(((0,0),(10,0),(10,10),(0,10),(0,6),(8,6),(8,5),(0,5)))

borderbuffered = border.buffer(-1)

clean1 = rdp(borderbuffered[0].exterior.coords[:-1],epsilon=0.5)
cleanpolygon1 = Polygon(clean1)
cleanpolygon1 #in loop, how to create these polygons depend on the length of the multipolygon? 
clean1


coords=(clean1[1],clean1[2])
(x1,y1),(x2,y2) = coords 
dx = x2 - x1
dy = y2 - y1
linelen = math.hypot(dx, dy)

x3 = x2 + dx/linelen * 0.1
y3 = y2 + dy/linelen * 0.1 

print(x3,y3)
