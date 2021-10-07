#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 23 16:08:44 2021

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
from shapely.ops import transform
import re
from rdp import rdp
import matplotlib.pyplot as plt
import copy
from collections import defaultdict
import os


border = Polygon(((0,0),(10,0),(10,10),(0,10),(0,6),(8,6),(8,5),(0,5)))
borderbuffered = border.buffer(-1)

TotalEdges = dict()

def addEdge(v, w):
    try:
        TotalEdges[v].add(w)
    except:
        TotalEdges[v]=set()
        print("create ",v)
        TotalEdges[v].add(w)
    try:
        TotalEdges[w].add(v)
    except:
        TotalEdges[w]=set()
        print("create ",w)
        TotalEdges[w].add(v)
        print("add ",v," to ",w)
        TotalEdges[v].add(w)
        print("add ",w," to ",v)
        print(TotalEdges)
        
borderCoord = border.exterior.coords[:-1]
BorderReflexPoints=[]
for i in range(len(borderCoord)):
    TotalEdges[borderCoord[i]] = set()
    addEdge(borderCoord[i],borderCoord[(i+1)%len(borderCoord)])
TotalEdges