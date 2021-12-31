import sys
import matplotlib.pyplot as plt
import itertools
from networkx.convert_matrix import _generate_weighted_edges
from networkx.generators.line import line_graph
from utl import *
import numpy as np
from rdp import rdp
import numpy.linalg as LA
import copy
from collections import defaultdict
import networkx as nx
from networkx.classes.function import path_weight
from shapely.geometry import *


TotalEdges = dict()
bufferpolygons = []
vertices=set()
BorderSet = set()
Graph = nx.MultiGraph()
radius=0.5

border = Polygon(((0, 3), (0, 7), (1, 7), (1, 8), (0, 8), (0, 14), (3, 14), (3, 11), (4, 11), (4, 14), (7, 14), (7, 8),
                  (4, 8), (4, 10), (3, 10), (3, 8), (2, 8), (2, 7), (3, 7), (3, 5), (5, 5), (5, 3), (6, 3), (6, 0),
                  (3, 0), (3, 3), (4, 3), (4, 4), (2, 4), (2, 3)))

borderbuffered = border.buffer(-radius)

def findpath(border):
    path=[]
    border=Polygon(border)
    borderbuffered = border.buffer(-radius)

    def findbisector(before):
        r=[]
        if True:
            for i in range(0,len(before)-1):
                u=[a-b for (a, b) in zip(before[i],before[i - 1])]
                v=[a-b for (a, b) in zip(before[i],before[i + 1]%len(before))] 
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

'''
G = nx.DiGraph()
G.add_edge((1,1),(2,2),parent=BufferedPolygon(kind="visible Polygon",
                        epsilon=0.1, radius=radius),weight=3)
print(G.edges(data=True))
'''


allpoints=set()

index=0
for i in borderbuffered[:]:

    a=BufferedPolygon(kind="visible Polygon",
                        item=i.exterior.coords[:-1],
                        epsilon=0.1, radius=radius, border=border)
    xBottom, yBottom = zip(*a.clean)
    plt.gca().plot(xBottom, yBottom, color="black", linestyle='--', linewidth=1.0) 
    bufferpolygons.append(a)
    coordY,_ = a.findbisector()
    a.addoutpoints(coordY)

    for i,coord in enumerate(coordY):
        vertices.add(Node(a,coord,index+i,"polygon"))
    
    for coord in a.clean:
        allpoints.add(tuple(coord))


    length=len(coordY)
    a.addrange([index,index+length])
    index=index+length

points=set()
bitangent(borderCoord=copy.deepcopy(border.exterior.coords[:-1]),border=border,points=points)
for i,coord in enumerate(points):
    vertices.add(Node(np.NaN,coord,index+i,"border")) #if parent=np.NaN -> root

for coord in border.exterior.coords[:-1]:
    allpoints.add(tuple(coord))

def fullnetwork():
    G=nx.MultiGraph()
    adjmatrix = [ [0]*len(allpoints) for i in range(len(allpoints))]
    for rowele in allpoints:
    #print(f"{i}: ({round(rowele.location[0],2)}, {round(rowele.location[1],2)})")
    #plt.text(rowele.location[0], rowele.location[1]+0.2, f'{i}',fontsize=12,fontname="Times New Roman")

        for columnele in allpoints:
            if (rowele==columnele):
                continue

            l=LineString([columnele,rowele])
    
            if l.within(border):
                G.add_edge(rowele,columnele,weight=l.length)

    return G

Gfull = fullnetwork()

def validnetwork():
    G = nx.DiGraph()
    adjmatrix = [ [0]*len(vertices) for i in range(len(vertices))]
    for i, rowele in enumerate(vertices):
        #print(f"{i}: ({round(rowele.location[0],2)}, {round(rowele.location[1],2)})")
        #plt.text(rowele.location[0], rowele.location[1]+0.2, f'{i}',fontsize=12,fontname="Times New Roman")

        for j, columnele in enumerate(vertices):
            if (i==j):
                continue

            l=LineString([tuple(rowele.location),tuple(columnele.location)])
    
            if l.within(border):
                G.add_edge(rowele,columnele,weight=l.length)

Gvalid=validnetwork()

def compressednetwork(G):
    g1 = nx.DiGraph()

    for i in range(index):
        for j in range(index):
            if (i==j):
                continue
            weight,path = nx.bidirectional_dijkstra(G,i,j)
            g1.add_edge(i,j,path=path,weight=weight)
    return g1

Gcompressed = compressednetwork(Gvalid)

path = nx.approximation.greedy_tsp(g1)
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


