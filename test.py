import sys

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

#print(os.path.dirname(os.path.abspath(__file__)))
#sys.exit()

print("here")
TotalEdges = dict()
BorderSet = set()

border = Polygon(((0, 3), (0, 7), (1, 7), (1, 8), (0, 8), (0, 14), (3, 14), (3, 11), (4, 11), (4, 14), (7, 14), (7, 8),
                  (4, 8), (4, 10), (3, 10), (3, 8), (2, 8), (2, 7), (3, 7), (3, 5), (5, 5), (5, 3), (6, 3), (6, 0),
                  (3, 0), (3, 3), (4, 3), (4, 4), (2, 4), (2, 3)))

borderbuffered = border.buffer(-1)


borderFOUR = Polygon(((0, 3), (0, 7), (1, 7), (1, 8), (0, 8), (0, 14), (3, 14), (3, 11), (4, 11), (4, 14), (7, 14),
                      (7, 8), (4, 8), (4, 10), (3, 10), (3, 8), (2, 8), (2, 7), (3, 7), (3, 5), (5, 5), (5, 3), (6, 3),
                      (6, 0), (3, 0), (3, 3), (4, 3), (4, 4), (2, 4), (2, 3)))
borderFOURbuffered = borderFOUR.buffer(-1)

# coord1 = copy.deepcopy(rdp(borderbuffered[0].exterior.coords[:],epsilon=0.5))
# coord2 = copy.deepcopy(rdp(borderbuffered[1].exterior.coords[:],epsilon=0.5))
coord3 = copy.deepcopy(borderFOUR.exterior.coords[:])

def graphnumberedmap():
    plt.figure(figsize=(3.5, 6), dpi=100)

    # xBottom, yBottom = zip(*coord1)
    # xTop, yTop = zip(*coord2)
    xBorder, yBorder = zip(*coord3)

    # for i_x, i_y in coord1:
    # plt.text(i_x, i_y, '({:0.2f}, {:0.2f})'.format(i_x, i_y))

    # for i_x, i_y in coord2:
    # plt.text(i_x, i_y, '({:0.2f}, {:0.2f})'.format(i_x, i_y)) #format the numbers

    # for i_x, i_y in coord3:
    # plt.text(i_x, i_y, '({:0.2f}, {:0.2f})'.format(i_x, i_y)) #format the numbers

    # plt.gca().plot(xBottom,yBottom)
    # plt.gca().plot(xTop,yTop)
    plt.gca().plot(xBorder, yBorder)
    # plt.title("after cleaning")
    x = 1
    for i in borderFOURbuffered:
        coordX = copy.deepcopy(rdp(i.exterior.coords[:], epsilon=0.1))  # 6655 verticies
        Px, Py = zip(*coordX)
        plt.gca().plot(Px, Py)
        plt.text(i.representative_point().coords[0][0], i.representative_point().coords[0][1], '({:0.2f}, {:0.2f})'.format(i.representative_point().coords[0][0],i.representative_point().coords[0][1]))
        x += 1

    plt.show()



def addEdge(v, w):
    try:
        TotalEdges[v].add(w)
    except:
        TotalEdges[v] = set()
        TotalEdges[v].add(w)
    try:
        TotalEdges[w].add(v)
    except:
        TotalEdges[w] = set()
        TotalEdges[w].add(v)
        TotalEdges[v].add(w)


totalPath = []

# create graph
Vertices = []

borderCoord = border.exterior.coords[:-1]

BorderReflexPoints = []
for i in range(len(borderCoord)):
    addEdge(borderCoord[i], borderCoord[(i + 1) % len(borderCoord)])
    BorderSet.add(borderCoord[i])




# find reflex point
visiblePolygon = []
Border = []
EdgeSet = dict()
visited = []






BorderReflexPoints = []
for i in range(len(borderCoord)):
    TotalEdges[borderCoord[i]] = set()
    if isConvex(borderCoord[i - 1], borderCoord[i], borderCoord[((i + 1) % (len(borderCoord)))]):
        BorderReflexPoints.append(borderCoord[i])

for i in range(len(borderbuffered)):
    BufferedPolygon(kind="visible Polygon",
                    item=borderbuffered[i].exterior.coords[:-1],
                    epsilon=0.5,
                    parent=visiblePolygon)

#for i in BorderReflexPoints:
  #  for j in BorderReflexPoints:
       # if i!=j:
            #if not (crossBorder(borderCoord,i, j)):


def visit(startPolygon: BufferedPolygon):
    startPolygon.peeked = True
    existDirect = False
    for otherPolygons in visiblePolygon:
        if otherPolygons.ID != startPolygon.ID:
            for start in startPolygon.reflexPoints:
                for others in otherPolygons.reflexPoints:
                    if not (crossBorder(borderCoord, start, others)):
                        existDirect = True
                        e = Edge(polygonStart=startPolygon, polygonDest=otherPolygons, path=[start, others])
                        if startPolygon.dontHave(e):
                            startPolygon.edge.append(e)
                        if otherPolygons.dontHave(e):
                            # print("add edge to other polygon")
                            otherPolygons.edge.append(e)
                        try:
                            EdgeSet[tuple(start)].append(e)
                        except:
                            EdgeSet[tuple(start)] = []

    for PointB in BorderReflexPoints:
        for start in startPolygon.reflexPoints:
            if not penetrateBorder(border,start, PointB):
                e = Edge(polygonStart=startPolygon, polygonDest=border, path=[start, PointB])
                if tuple(start) in EdgeSet.keys():
                    EdgeSet[tuple(start)].append(e)
                else:
                    EdgeSet[tuple(start)]= []
                    EdgeSet[tuple(start)].append(e)
                TotalEdges[tuple(start)] = set()
                TotalEdges[PointB] = set()
                TotalEdges[tuple(start)].add(PointB)
                TotalEdges[PointB].add(tuple(start))

    if existDirect:
        sortedEdges = sorted(startPolygon.edge)
        for edge in sortedEdges:
            if not edge.used and not edge.polygonDest.peeked:
                # if edge.polygonDest!=startPolygon: #if there exist an edge with the current polygon as dest, most probably its already peeked or visite
                polygonNext = edge.polygonDest
                if not polygonNext.peeked:
                    edge.used = True
                    totalPath.append(edge)
                    visit(polygonNext)

    else:
        # use border
        for start in startPolygon.reflexPoints:
            queue = []
            queue.append(start)

            # while queue:
            # s = queue.pop(0)

            # for i in self.graph[s]:
            # if i not in BorderSet:
            # queue.append(i)



#visit(visiblePolygon[0])


    
def graphinnermap():
    for i in borderbuffered:
        store=[]
        reduce(poly=i, epsilon=0.1, buffer=-0.3, areamin=0.5, store=store)
        for c in store:
            xBottom, yBottom = zip(*c)
            plt.gca().plot(xBottom, yBottom, color="black", linestyle='--', linewidth=1.0)  # marker=''

    coord1 = copy.deepcopy(border.exterior.coords[:])
    xBottom, yBottom = zip(*coord1)
    plt.gca().plot(xBottom, yBottom, color="black", linewidth=1.0)

    # plt.title("after cleaning")
    plt.savefig('C:/Users/space/Documents/comp400 files/contour.png')
    plt.show()


def graphbitangentmap():
    plt.figure(figsize=(3.5, 6), dpi=200)

    coord3 = copy.deepcopy(borderFOUR.exterior.coords[:])
    xBorder, yBorder = zip(*coord3)
    plt.gca().plot(xBorder, yBorder,color='black',linewidth=1.0)

    #line = LineString(borderFOUR.exterior.coords[:])


    for i in borderbuffered:
        coord1 = rdp(list(copy.deepcopy(i.exterior.coords[:])), epsilon=0.1)
        xBottom, yBottom = zip(*coord1)
        plt.gca().plot(xBottom, yBottom, color="black", linestyle='--', linewidth=1.0)  # marker=''
        b2=i.buffer(-1)

    s1=[]
    s2=[]
    bitangent(borderCoord=borderCoord,border=border,store1=s1,store2=s2)
    for i in s1:
        Px, Py = zip(*i) #zip is useful
        plt.gca().plot(Px, Py, color='black',linestyle='-.',marker='o',linewidth=1.0)
    for i in s2:
        Px, Py = zip(*i) #zip is useful
        plt.gca().plot(Px, Py, color='black',linestyle=':',marker='o',linewidth=2.0) 

    plt.savefig('C:/Users/space/Documents/comp400 files/map2.png')
    plt.show() 




'''
for polygons in visiblePolygon:
    #find the path from the last polygon to this one
    if polygons.peeked == False:
        print("visit", polygons)
        visit(polygons)
print("done!")
 
print(totalPath)
'''
