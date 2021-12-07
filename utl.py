
#from shapely.geometry import Point, Polygon,LineString
import matplotlib.pyplot as plt
import itertools
from utl import *
# from shapely.geometry import LineString,mapping,shape,MultiLineString,MultiPolygon,Polygon,MultiPoint,LinearRing
# import json
# import os.path
# from svgpath2mpl import parse_path
import numpy as np
# import matplotlib.pyplot as plt
# from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, Close
# from matplotlib.path import Path
# import matplotlib.patches as patches
# from shapely.ops import transform
# from shapely.affinity import translate
# import re
# import textwrap
from rdp import rdp
# import matplotlib.pyplot as plt
import math
import numpy.linalg as LA
import copy
from collections import defaultdict
# import itertools
import os
from shapely.geometry import *



def reduce(poly:Polygon,epsilon,buffer=1,areamin=0,store=[]):
    coord = rdp(list(copy.deepcopy(poly.exterior.coords[:])), epsilon)
    if len(coord)>=2:
        #print("area to graph",poly.area)
        store.append(coord)
        b2=poly.buffer(buffer)
        if b2.area>=areamin:
            reduce(b2,epsilon,buffer,areamin,store) 
    return



def bitangent(borderCoord,border:Polygon,store1=[],store2=[]): #why use this for convex: avoid corner cases
    line = LineString(border.exterior.coords[:])
    for (x, y) in itertools.combinations(borderCoord, 2):
        if not crossBorder(borderCoord,x, y) and not penetrateBorder(border,y, x) and not penetrateBorder(border,x, y):
            if Point(newcoord((x, y), -0.1)).within(border) and Point((x[0] + y[0]) / 2, (x[1] + y[1]) / 2).within(border):
                coordX = [x, y]
                #Px, Py = zip(*coordX) #zip is useful
                store1.append(coordX)
            
            elif line.contains(Point(newcoord((x, y), -0.1))) and line.contains(
                Point((x[0] + y[0]) / 2, (x[1] + y[1]) / 2)):
                coordX = [x, y]
                store2.append(coordX)
                
    

def onSegment(p, q, r):
    if ((q[0] <= max(p[0], r[0])) and (q[0] >= min(p[0], r[0])) and
            (q[1] <= max(p[1], r[1])) and (q[1] >= min(p[1], r[1]))):
        return True
    return False


def orientation(p, q, r):
    val = (float(q[1] - p[1]) * float(r[0] - q[0])) - (float(q[0] - p[0]) * float(r[1] - q[1]))
    if (val > 0):
        return 1
    elif (val < 0):
        return 2
    else:
        return 0


def doCross(p1, q1, p2, q2):  # touching or parallel do not cross, only strictly intersect count as cross
    o1 = orientation(p1, q1, p2)
    o2 = orientation(p1, q1, q2)
    o3 = orientation(p2, q2, p1)
    o4 = orientation(p2, q2, q1)

    if ((o1 == 0) and onSegment(p1, p2, q1)):
        return False
    if ((o2 == 0) and onSegment(p1, q2, q1)):
        return False
    if ((o3 == 0) and onSegment(p2, p1, q2)):
        return False
    if ((o4 == 0) and onSegment(p2, q1, q2)):
        return False
    if ((o1 != o2) and (o3 != o4)):
        return True
    return False


def crossBorder(borderCoord, p, q):
    for i in range(len(borderCoord)):
        # print(borderCoord[i-1],borderCoord[i]) #counter clockwise
        val = doCross(borderCoord[i - 1], borderCoord[i], p, q)
        if val:
            return True


def newcoord(coords, dist):  # coords: (A,B) -- from A to B. Extend after B.      #啊啊啊python真的太难受了！！x1是什么，是什么type啊！！
    (x1, y1), (x2, y2) = coords
    dx = x2 - x1
    dy = y2 - y1
    linelen = math.hypot(dx, dy)

    x3 = x2 + dx / linelen * dist
    y3 = y2 + dy / linelen * dist
    return x3, y3


def penetrateBorder(border,Pstart, Pend):
    return not Point(newcoord((Pstart, Pend), 0.1)).within(
        border)  # within means included within. so penetrate means outside


def isConvex(p1, p2, p3):  # return true if convex or straight
    a = np.array(p1)
    b = np.array(p2)
    c = np.array(p3)
    ab = a - b
    cb = c - b
    cross = np.cross(cb,
                     ab)  # negative degree -> convex 凸                                                         #https://stackoverflow.com/questions/20252845/best-algorithm-for-detecting-interior-and-exterior-angles-of-an-arbitrary-shape
    dot = np.inner(cb, ab)
    radtan = np.arctan2(cross, dot)
    degtan = np.rad2deg(radtan)
    return degtan <= 0

    
def pathDist(path=[]):
    dist = 0
    if (len(path) == 0):
        raise ValueError("empty path")
    else:
        for i in range(len(path) - 1):
            a = np.array(path[i])
            b = np.array(path[i + 1])
            dist += np.linalg.norm(a - b)
            # print (path[i],path[i+1],dist)
    return dist



class Edge:
    def __init__(self, polygonStart, polygonDest, path=[]):
        try:
            self.polygonStart = polygonStart
            self.polygonDest: BufferedPolygon = polygonDest
            self.used = False
            self.path = path
            self.weight = pathDist(path)
        except:
            print("exception catched: empty path")

    def getWeight(self):
        return self.weight

    def appendPath(self, p):
        self.path.append(p)

    def isSameAs(self, other):
        locationSame1: bool = self.polygonStart == other.polygonStart
        locationSame2: bool = self.polygonDest == other.polygonDest
        weightSame: bool = self.weight == other.weight
        pathSame: bool = self.path == other.path
        return (pathSame and weightSame and locationSame1 and locationSame2)

    def __eq__(self, other):
        weightSame = self.weight == other.weight
        return (weightSame)

    def __ne__(self, other):
        return (self.weight != other.weight)

    def __lt__(self, other):
        return (self.weight < other.weight)

    def __le__(self, other):
        return (self.weight <= other.weight)

    def __gt__(self, other):
        return (self.weight > other.weight)

    def __ge__(self, other):
        return (self.weight >= other.weight)

    def __repr__(self):
        return "Edge from {} to {} with path {} of length {} \n".format(
            self.polygonStart, self.polygonDest, self.path, self.weight) + os.linesep

    def __str__(self):
        return "Edge from {} to {} with path {} of length {} \n".format(
            self.polygonStart, self.polygonDest, self.path, self.weight) + os.linesep


# container for buffered polygon storing all the information needed
class BufferedPolygon:

    def __init__(self, kind="", item=[(0, 1), (1, 1), (1, 0)], epsilon=0, parent=[]):
        self.raw = item
        self.clean = rdp(item, epsilon=epsilon)
        self.reflexPoints = self.findreflex(self.clean)
        self.name = kind
        self.parent = parent
        parent.append(self)
        self.ID = Polygon(item).representative_point().coords[:]
        self.edge = []
        self.visited = False
        self.peeked = False

    def dontHave(self, e: Edge):
        if len(self.edge) == 0:
            return True
        else:
            for thisE in self.edge:
                if thisE.isSameAs(e):
                    return False
            return True

    def add(self, reflex):
        self.reflexPoints.append(reflex)

    def findreflex(self, before):
        reflexPoints = []
        for i in range(len(before)):
            if isConvex(before[i - 1], before[i], before[((i + 1) % (len(before)))]):
                reflexPoints.append(
                    before[i])  # this "self" thing.. I want java when you can just refer to variables by their name :(
        return reflexPoints  # 啊我写python真的就是连猜带蒙啊。。

    def __repr__(self):
        return "{} at {} \n".format(self.name, self.ID) + os.linesep

    def __str__(self):
        return "{} at {} \n".format(self.name, self.ID)