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
radius=1

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

   # plt.show()

'''
G=nx.MultiGraph()
#G.add_edge(1,0,weight=10)
G.add_edge(2,1,weight=12)
G.add_edge(1,0,weight=10)
G.add_edge(3,0,weight=3)
m = nx.to_numpy_matrix(G)

print(m)
print(path_weight(G, [2,1,0],weight="weight"))
sys.exit()
'''

G = nx.DiGraph()
vertices=[]
bufferpolygons = []
index=0

for i in borderbuffered[:]:
    a=BufferedPolygon(kind="visible Polygon",
                        item=i.exterior.coords[:-1],
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

points=set()
bitangent(borderCoord=copy.deepcopy(border.exterior.coords[:-1]),border=border,points=points)
for i,coord in enumerate(points):
    vertices.append(Node(a,coord,index+i,"border"))



adjmatrix = [ [0]*len(vertices) for i in range(len(vertices))]
for i, rowele in enumerate(vertices):
    #print(f"{i}: ({round(rowele.location[0],2)}, {round(rowele.location[1],2)})")
    #plt.text(rowele.location[0], rowele.location[1]+0.2, f'{i}',fontsize=15,fontname="Times New Roman")

    for j, columnele in enumerate(vertices):
        if (i==j):
            continue

        l=LineString([tuple(rowele.location),tuple(columnele.location)])
 
        if l.within(border):
            adjmatrix[i][j]=l.length

#print('\n'.join([' '.join(['{:0.2f}'.format(item) for item in row]) 
     # for row in adjmatrix]))

#plt.show()

G = nx.DiGraph(np.array(adjmatrix))
#compressed path
g1 = nx.DiGraph()
gtest=nx.DiGraph()
gtest.add_node(0)
for i in range(index):
    for j in range(index):
        if (i==j):
            continue
        weight,path = nx.bidirectional_dijkstra(G,i,j)
      #  g1.add_node(1)
      #  g1.add_node(2)
      #  print(g1.nodes)
        g1.add_edge(i,j,path=path,weight=weight)
        gtest.add_edge(0,i+1,path=[0,i+1],weight=0)
        gtest.add_edge(i+1,j+1,path=path,weight=weight)

#print(path_weight(gtest,[0,3,2],weight="weight"))
g1.add_edge(6,2,weight=0)
g1.add_edge(5,6,weight=0)
tsp = nx.approximation.traveling_salesman_problem
path = tsp(g1)
print(path)


sys.exit()
tsp = nx.approximation.traveling_salesman_problem
path = tsp(gtest)
print(path)
sys.exit()




g1.add_node(6)
g1.add_edge(6,0,weight=0)
g1.add_edge(6,1,weight=0)
g1.add_edge(6,2,weight=0)
g1.add_edge(6,3,weight=0)
g1.add_edge(6,4,weight=0)
g1.add_edge(6,5,weight=0)
tsp = nx.approximation.traveling_salesman_problem
path = tsp(g1)
print(printfullpath(g1,path))
sys.exit()
'''
g1.add_edge(6,0,weight=0)
weight,path = nx.bidirectional_dijkstra(g1,6,0) #0,6 dont exist
print(g1.get_edge_data(*path,"path")['path'])
SA_tsp = nx.approximation.greedy_tsp(g1, source=)
#method = lambda g1, wt: SA_tsp(g1, "greedy", weight=wt, temp=500)
#tsp = nx.approximation.traveling_salesman_problem
tsp(g1, nodes=range(index))
path = tsp(g1, method=method)
print(path)
p=g1.get_edge_data(*path,"path")['path']
print(path_weight(g1, path,weight="weight"))
pweight= path_weight(g1, [2,1,0,3,4,5],weight="weight")
'''


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

path,weight = partialpath(g1)
print(path,weight)
sys.exit()



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
    Graph.add_node(BufferedPolygon(kind="visible Polygon",
                    item=borderbuffered[i].exterior.coords[:-1],
                    epsilon=0.5,
                    parent=visiblePolygon))

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

    plt.savefig('C:/Users/space/Documents/comp400 files/map2.png')
    plt.show() 


#graphinnermap()

'''
for polygons in visiblePolygon:
    #find the path from the last polygon to this one
    if polygons.peeked == False:
        print("visit", polygons)
        visit(polygons)
print("done!")
 
print(totalPath)
'''
