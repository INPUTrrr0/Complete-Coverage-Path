import enum
import sys

#from shapely.geometry import Point, Polygon,LineString
import matplotlib.pyplot as plt
import itertools
from networkx.convert_matrix import _generate_weighted_edges

from networkx.generators.line import line_graph
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
import networkx as nx
from networkx.classes.function import path_weight
from shapely.geometry import *

#print(os.path.dirname(os.path.abspath(__file__)))
#sys.exit()


TotalEdges = dict()
BorderSet = set()
Graph = nx.MultiGraph()
radius=0.5

border = Polygon(((0, 3), (0, 7), (1, 7), (1, 8), (0, 8), (0, 14), (3, 14), (3, 11), (4, 11), (4, 14), (7, 14), (7, 8),
                  (4, 8), (4, 10), (3, 10), (3, 8), (2, 8), (2, 7), (3, 7), (3, 5), (5, 5), (5, 3), (6, 3), (6, 0),
                  (3, 0), (3, 3), (4, 3), (4, 4), (2, 4), (2, 3)))

borderbuffered = border.buffer(-radius)

'''
a=[1,2]
b=[0,0]
l=LineString([tuple(a),tuple(b)])

sys.exit()
'''

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

def findbisector(before):
    r=[]
    if True:
        for i in range(1,len(before)-1):
            u=[a-b for (a, b) in zip(before[i],before[i - 1])]
            v=[a-b for (a, b) in zip(before[i],before[i + 1] )] 
            Mu = LA.norm(u)
            Mv = LA.norm(v)
            u = [number / Mu for number in u]
            v = [number / Mv for number in v]
            b= [a+b for (a, b) in zip(u,v)]
            #b=[i.tolist() for i in b]
            new = [x+radius*y for (x, y) in zip(before[i],b)]
            newP=Point(new)

            inout = border.contains(newP)
            if (inout):
                r.append(b)
    return r


def validPoints():
    plt.figure(figsize=(3.5, 6), dpi=100)
    coord1 = copy.deepcopy(border.exterior.coords[:])
    xBottom, yBottom = zip(*coord1)
    plt.gca().plot(xBottom, yBottom, color="black", linewidth=1.0)
    # plt.title("after cleaning")

    for i in borderbuffered[:]:
        a=BufferedPolygon(kind="visible Polygon",
                        item=i.exterior.coords[:-1],
                        epsilon=0.1, radius=radius, border=border)
        xBottom, yBottom = zip(*a.clean)
        plt.gca().plot(xBottom, yBottom, color="black", linestyle='--', linewidth=1.0)  # marker=''
        coordX = a.reflexPoints
        coordY,_ = a.findbisector()
        for r in coordY:
            plt.plot(r[0],r[1],'kv',markersize=10) 
        s1=[]
        s2=[]
        points=set()
        bitangent(borderCoord=coord1[:-1],border=border,store1=s1,store2=s2,points=points)
        for p in points:
            plt.plot(p[0],p[1],'k.',markersize=10) 
    plt.show()


G = nx.DiGraph()
vertices=[]
bufferpolygons = []
index=0


xBorder, yBorder = zip(*border.exterior.coords)
plt.gca().plot(xBorder, yBorder,color="blue")

store=[]
for i in borderbuffered[:]:
    reduce(border,epsilon=0.1,buffer=-radius,areamin=1,store=store)
store=[ele for ind, ele in enumerate(store) if ele not in store[:ind]]
#print(len(store))

'''
for i in borderbuffered[:]:
    xBottom, yBottom = zip(*i)
    plt.gca().plot(xBottom, yBottom, color="black", linestyle='--', linewidth=1.0) 
    a=BufferedPolygon(kind="visible Polygon",
                        item=i,
                        epsilon=0.1, radius=radius, border=border)
    xBottom, yBottom = zip(*a.clean)
    bufferpolygons.append(a)
    coordY,_ = a.findbisector()
    a.addoutpoints(coordY)

    for i,coord in enumerate(coordY):
        vertices.append(Node(a,coord,index+i,"polygon"))

    length=len(coordY)
    a.addrange([index,index+length])
    index=index+length
'''

for i in borderbuffered[:]:
    #xBottom, yBottom = zip(*i)
    #plt.gca().plot(xBottom, yBottom, color="black", linestyle='--', linewidth=1.0) 
    a=BufferedPolygon(kind="visible Polygon",
                        item=i.exterior.coords,
                        epsilon=0.1, radius=radius, border=border)
    xBottom, yBottom = zip(*a.clean)
    plt.gca().plot(xBottom, yBottom, color="black", linestyle='--', linewidth=1.0) 
    bufferpolygons.append(a)
    coordY,_ = a.findbisector()
    a.addoutpoints(coordY)

    for i,coord in enumerate(a.clean):
        vertices.append(Node(a,coord,index+i,"polygon"))

    length=len(a.clean)
    a.addrange([index,index+length])
    index=index+length

points=set()
#bitangent(borderCoord=copy.deepcopy(border.exterior.coords[:-1]),border=border,points=points)
for i,coord in enumerate(border.exterior.coords[:-1]):
    vertices.append(Node(a,coord,index+i,"border"))



adjmatrix = [ [0]*len(vertices) for i in range(len(vertices))]
for i, rowele in enumerate(vertices):
    print(f"{i}: ({round(rowele.location[0],2)}, {round(rowele.location[1],2)})")
    plt.text(rowele.location[0], rowele.location[1]+0.2, f'{i}',fontsize=12,fontname="Times New Roman")

    for j, columnele in enumerate(vertices):
        if (i==j):
            continue

        l=LineString([tuple(rowele.location),tuple(columnele.location)])
 
        if l.within(border):
            adjmatrix[i][j]=l.length

#print('\n'.join([' '.join(['{:0.2f}'.format(item) for item in row])])

     # for row in adjmatrix]))


G = nx.DiGraph(np.array(adjmatrix))

g1 = nx.DiGraph()


for i in range(len(vertices)):
    for j in range(len(vertices)):
        if (i==j):
            continue
        weight,path = nx.bidirectional_dijkstra(G,i,j)
        g1.add_edge(i,j,path=path,weight=weight)

path = nx.approximation.greedy_tsp(g1,source=2)
pweight= path_weight(g1, path,weight="weight")
print(path,pweight)
plt.show()
sys.exit()


#compressed path without border
g1 = nx.DiGraph()


for i in range(index):
    for j in range(index):
        if (i==j):
            continue
        weight,path = nx.bidirectional_dijkstra(G,i,j)
        g1.add_edge(i,j,path=path,weight=weight)



def fullpath(g1):
    shortestpath=[]
    weight = float('inf') 
    for i in range(index):
        path = nx.approximation.greedy_tsp(g1,source=i)
        pweight= path_weight(g1, path,weight="weight")
        if pweight<weight:
            weight=pweight
            shortestpath=path
    return printfullpath(g1,shortestpath),weight


def partialpath(g1):
    shortestpath=[]
    weight = float('inf') 
    adjtemp = nx.to_numpy_matrix(g1)
    for _,ele in enumerate(bufferpolygons):
        for i in range(*ele.range):
            for j in [x for x in range(index) if x not in range(*ele.range)]:
                gtemp=nx.DiGraph(adjtemp)
                gtemp.add_edge(j,i,weight=0) #from dest to start, weight=0
                path = nx.approximation.greedy_tsp(gtemp, source=i)
                pweight= path_weight(gtemp, path,weight="weight")
                if pweight<weight:
                    weight=pweight
                    shortestpath=path
    return printfullpath(g1,shortestpath[:-1]),weight

path,weight = fullpath(g1)
print(path,weight)
sys.exit()


    
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
    #plt.savefig('C:/Users/space/Documents/comp400 files/contour.png')
    plt.show()


def graphbitangentmap():
    plt.figure(figsize=(3.5, 6), dpi=200)

    coord3 = copy.deepcopy(border.exterior.coords[:])
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

    #plt.savefig('C:/Users/space/Documents/comp400 files/map2.png')
    plt.show() 


