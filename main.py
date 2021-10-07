'''
Created on Aug. 23, 2021

@author: theone
'''

from shapely.geometry import Point,Polygon
#from shapely.geometry import LineString,mapping,shape,MultiLineString,MultiPolygon,Polygon,MultiPoint,LinearRing
#import json
#import os.path
#from svgpath2mpl import parse_path
import numpy as np 
#import matplotlib.pyplot as plt
#from svg.path import Path, Line, Arc, CubicBezier, QuadraticBezier, Close
#from matplotlib.path import Path
#import matplotlib.patches as patches
#from shapely.ops import transform
#from shapely.affinity import translate
#import re
#import textwrap
from rdp import rdp
#import matplotlib.pyplot as plt
import math
import numpy.linalg as LA
import copy
from collections import defaultdict
#import itertools
import os


TotalEdges = dict()
BorderSet = set()

border = Polygon(((0,3),(0,7),(1,7),(1,8),(0,8),(0,14),(3,14),(3,11),(4,11),(4,14),(7,14),(7,8),(4,8),(4,10),(3,10),(3,8),(2,8),(2,7),(3,7),(3,5),(5,5),(5,3),(6,3),(6,0),(3,0),(3,3),(4,3),(4,4),(2,4),(2,3)))
borderbuffered = border.buffer(-1)


def addEdge(v, w):
    try:
        TotalEdges[v].add(w)
    except:
        TotalEdges[v]=set()
        TotalEdges[v].add(w)
    try:
        TotalEdges[w].add(v)
    except:
        TotalEdges[w]=set()
        TotalEdges[w].add(v)
        TotalEdges[v].add(w)


totalPath=[]

#create graph
Vertices = []
    
borderCoord = border.exterior.coords[:-1]

BorderReflexPoints=[]
for i in range(len(borderCoord)):
    addEdge(borderCoord[i],borderCoord[(i+1)%len(borderCoord)])
    BorderSet.add(borderCoord[i])
    
    
def onSegment(p, q, r):
    if ( (q[0] <= max(p[0], r[0])) and (q[0] >= min(p[0], r[0])) and 
           (q[1] <= max(p[1], r[1])) and (q[1] >= min(p[1], r[1]))):
        return True
    return False

def orientation(p,q,r):
    val = (float(q[1]-p[1])*float(r[0]-q[0]))-(float(q[0]-p[0])*float(r[1]-q[1]))
    if (val>0):
        return 1
    elif (val<0):
        return 2
    else:
        return 0
    
def doCross(p1,q1,p2,q2): # touching or parallel do not cross, only strictly intersect count as cross 
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

def crossBorder(p,q):
    for i in range(len(borderCoord)):
        #print(borderCoord[i-1],borderCoord[i]) #counter clockwise
        val = doCross(borderCoord[i-1],borderCoord[i],p,q)
        if val:
            return True
        
def newcoord(coords, dist): #coords: (A,B) -- from A to B. Extend after B.      #啊啊啊python真的太难受了！！x1是什么，是什么type啊！！
    (x1,y1),(x2,y2) = coords 
    dx = x2 - x1
    dy = y2 - y1
    linelen = math.hypot(dx, dy)

    x3 = x2 + dx/linelen * dist
    y3 = y2 + dy/linelen * dist    
    return x3, y3

def penetrateBorder(Pstart,Pend): 
    return not Point(newcoord((Pstart,Pend),0.1)).within(border) #within means included within. so penetrate means outside


def isConvex(p1,p2,p3): #return true if convex or straight 
    a = np.array(p1)
    b = np.array(p2)
    c = np.array(p3) 
    ab=a-b
    cb=c-b
    cross = np.cross(cb,ab) #negative degree -> convex 凸                                                         #https://stackoverflow.com/questions/20252845/best-algorithm-for-detecting-interior-and-exterior-angles-of-an-arbitrary-shape
    dot = np.inner(cb,ab)
    radtan = np.arctan2(cross,dot)
    degtan = np.rad2deg(radtan)
    return degtan<=0

#find reflex point
visiblePolygon = []
Border = []
EdgeSet=dict()
visited=[]

def pathDist(path=[]):
    dist=0
    if (len(path)==0):
        raise ValueError("empty path")
    else:
        for i in range(len(path)-1):
            a=np.array(path[i])
            b=np.array(path[i+1])
            dist += np.linalg.norm(a-b)
            #print (path[i],path[i+1],dist)
    return dist


class Edge:
    def __init__(self, polygonStart,polygonDest,path=[]):
        try:
            self.polygonStart = polygonStart
            self.polygonDest:BufferedPolygon = polygonDest
            self.used = False
            self.path = path
            self.weight = pathDist(path)
        except:
            print("exception catched: empty path")
    
    def getWeight(self):
        return self.weight
    
    def appendPath(self,p):
        self.path.append(p)
        
    def isSameAs(self,other):
        locationSame1: bool  = self.polygonStart == other.polygonStart
        locationSame2: bool  = self.polygonDest == other.polygonDest
        weightSame: bool =self.weight == other.weight
        pathSame: bool = self.path == other.path
        return (pathSame and weightSame and locationSame1 and locationSame2)
        
    def __eq__(self, other):
        weightSame=self.weight == other.weight
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
             self.polygonStart,self.polygonDest,self.path,self.weight)+os.linesep
    
    def __str__(self):
         return "Edge from {} to {} with path {} of length {} \n".format(
             self.polygonStart,self.polygonDest,self.path,self.weight)+os.linesep



#container for buffered polygon storing all the information needed
class BufferedPolygon:

    def __init__(self, kind="",item=[(0,1),(1,1),(1,0)],epsilon=0,parent=[]):
        self.raw=item
        self.clean=rdp(item,epsilon=epsilon)
        self.reflexPoints=self.findreflex(self.clean)
        self.name = kind
        self.parent = parent
        parent.append(self)
        self.ID = Polygon(item).representative_point().coords[:]
        self.edge=[]
        self.visited=False
        self.peeked=False
        
    def dontHave(self, e:Edge):
            if len(self.edge)==0:
                    return True
            else:
                for thisE in self.edge:
                    if thisE.isSameAs(e):
                        return False
                return True

        
    def add(self,reflex):
        self.reflexPoints.append(reflex)
        
    def findreflex(self,before):
        reflexPoints=[]
        for i in range(len(before)):
            if isConvex(before[i-1],before[i],before[((i+1)%(len(before)))]):
                reflexPoints.append(before[i])        #this "self" thing.. I want java when you can just refer to variables by their name :(
        return reflexPoints #啊我写python真的就是连猜带蒙啊。。
            
    def __repr__(self):
         return "{} at {} \n".format(self.name,self.ID)+os.linesep
    def __str__(self):
         return "{} at {} \n".format(self.name,self.ID)


BorderReflexPoints=[]
for i in range(len(borderCoord)):
    TotalEdges[borderCoord[i]] = set()
    if isConvex(borderCoord[i-1],borderCoord[i],borderCoord[((i+1)%(len(borderCoord)))]):
        BorderReflexPoints.append(borderCoord[i])   
        
        
for i in range(len(borderbuffered)):
    BufferedPolygon(kind="visible Polygon",
                    item=borderbuffered[i].exterior.coords[:-1],
                    epsilon=0.5,
                    parent=visiblePolygon)
    

def visit(startPolygon:BufferedPolygon): 
    startPolygon.peeked=True
    existDirect=False
    for otherPolygons in visiblePolygon:
        if otherPolygons.ID != startPolygon.ID:
            for start in startPolygon.reflexPoints:
                for others in otherPolygons.reflexPoints:
                    if not (crossBorder(start,others)):
                        existDirect=True
                        e=Edge(polygonStart=startPolygon,polygonDest=otherPolygons,path=[start,others])
                        if startPolygon.dontHave(e):
                                startPolygon.edge.append(e)     
                        if otherPolygons.dontHave(e):
                            #print("add edge to other polygon")
                            otherPolygons.edge.append(e)
                        try:
                            EdgeSet[tuple(start)].append(e)
                        except:
                            EdgeSet[tuple(start)]=[]
                            
                    
                        
            
    for PointB in BorderReflexPoints:
        for start in startPolygon.reflexPoints:
            if not penetrateBorder(start,PointB):
                e=Edge(polygonStart=startPolygon,polygonDest=border,path=[start,PointB])
                EdgeSet[tuple(start)].append(e)
                TotalEdges[tuple(start)]=set()
                TotalEdges[PointB]=set()
                TotalEdges[tuple(start)].add(PointB)
                TotalEdges[PointB].add(tuple(start))
    
    if existDirect:
        sortedEdges = sorted(startPolygon.edge)
        for i in range(len(sortedEdges)+1):
            edge = sortedEdges[i%len(sorted(startPolygon.edge))]
            if not edge.used:
                #if edge.polygonDest!=startPolygon: #if there exist an edge with the current polygon as dest, most probably its already peeked or visite
                    if not edge.polygonDest.peeked:
                        break
             
        polygonNext = edge.polygonDest
        if not polygonNext.peeked:
            edge.used=True
            totalPath.append(edge)
            visit(polygonNext)
    else:
        #use border 
        for start in startPolygon.reflexPoints:
            queue = []
            queue.append(start)
            
            #while queue:
                #s = queue.pop(0)

                #for i in self.graph[s]:
                    #if i not in BorderSet:
                       #queue.append(i)
        
        

visit(visiblePolygon[0])

for polygons in visiblePolygon:
        if polygons.peeked==False:
            visit(polygons)
            print("visit",polygons)
print("done!")
print(totalPath)
