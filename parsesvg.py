from svgpath2mpl import parse_path
import matplotlib.pyplot as plt
from rdp import rdp
import sys
from utl import *
import numpy as np
from shapely.geometry import Polygon, LineString
holes=[]

bordersvg= "m 2.12891,421.292 -0.90977,65.516 159.97986,-0.441 0.442,66.73 63.375,-1.809 0.156,-4.218 5.937,-0.313 0.313,5.625 18.281,-0.312 -3.281,129.687 9.375,0 0,13.125 16.875,0.625 0,7.5 52.187,0.625 0.313,-42.5 18.125,0.11 0.625,-21.36 7.5,-4.375 8.125,0 -1.25,55.625 34.258,0.735 0.261,-7.489 62.383,1.156 0.442,-14.402 15.687,0.195 -0.769,18.774 9.504,0.441 0.66,-28.285 -22.098,0.223 0,-6.629 20.992,0.441 1.18,-102.66 -16.875,-0.039 0.035,-7.316 7.395,0.261 1.632,-72.593 60.469,0 -0.469,25.625 46.875,0.625 0.313,-25.157 62.969,0 0,-17.187 8.566,-0.274 0,35.707 31.59,26.286 37.187,0.312 30.868,-21.98 1.32,-37.707 11.875,0 -0.625,11.875 50.625,-0.625 -1.875,67.187 17.344,0 -0.625,20.781 -14.219,0.157 -1.668,135.613 50.164,-0.223 0.555,-38.449 8.937,-0.937 -0.48,35.851 115.836,0.664 0,-19.668 12.814,-0.441 0,-43.309 -13.256,0.442 0,-27.403 -115.789,0.887 -1.09,40.398 -6.445,0.153 0.891,-83.266 -10,-0.156 0,-20.313 16.093,0.157 0.313,-5.625 21.562,0.312 0,23.438 65.625,-0.938 -1.875,-22.187 -44.062,0.937 -0.313,-7.187 45.313,0 -0.609,-23.946 -63.637,-0.996 0.219,22.43 -23.864,-0.551 3.516,-103.812 -109.063,1.25 -0.312,24.375 -10.313,0 2.704,-154.274 53.55,-0.262 -0.074,-15.207 12.813,0 -0.442,19.446 59.664,0 1.324,-129.93 49.278,0.328 -0.223,-13.035 -23.973,-0.664 -0.882,-18.891 -24.2,-0.222 -0.218,-22.981 24.527,-0.664 -0.223,13.039 42.207,0 0,-12.156 7.735,-0.219 0.219,-13.2581 9.945,-0.2219 0.219,-3.7539 28.066,0.4411 13.04,26.5158 29.6,12.152 0,21.875 -24.74,-0.66 -1.33,41.543 17.98,0.434 0.31,12.187 29.37,0.313 -0.31,-11.25 19.06,-1.25 0.32,-41.563 -25.94,0 -0.94,-20.312 31.56,-12.5 15.32,-23.438 6.32,0 -0.07,25.738 18.75,0.512 0,-7.5 11.25,0 -0.08,8.125 30.08,0 0,-6.32 11.14,0.07 1.36,-54.375 -55,0 -0.63,11.5629 -8.12,-0.3129 -0.63,-9.375 -8.12,0 -0.63,18.125 -5.62,0 -18.75,-28.125 -7.5,-3.75 2.39,-6.289 5.11,1.914 10.77,-22.6441 -17.02,-7.3559 -11.25,21.875 5,3.125 -1.88,6.25 -24.37,-8.75 8.6,-12.98 -8.09,-3.7821 -7.61,12.868 -13.89,-4.1207 -2.72,7.9019 14.4,5.3707 -22.57,12.8672 -16.873,23.75 -25,2.5 -1.875,-6.875 -10,0 2.5,-11.875 -132.5,-1.875 -0.281,20.4489 -11.309,-0.25 0,-3.9219 -19.227,-0.277 -0.441,15.027 20.551,0.223 0.223,-4.6398 9.214,0.2848 -1.039,171.406 -13.148,0.332 -0.11,-13.254 -54.8,-0.777 0.883,-41.985 -38.45,-0.882 0,-21.657 38.891,-0.883 0.293,-127.945 -22.813,0.4692 0.313,-43.5942 -15.625,-0.625 0.469,-9.68708 -21.875,-0.46875 0.156,10.46873 -16.563,0 -1.066,43.9571 -20.496,-0.2071 -2.672,128.0741 37.563,-0.129 -0.364,22.227 -39.41,0 -2.773,239.542 -9.063,-0.027 0,-15 -63.281,-0.625 -0.313,-26.875 -46.25,-0.625 0,27.813 -61.875,-0.313 -0.312,-10.625 -25,-0.937 0,17.187 -50,-1.25 -0.313,15.625 -7.5,0.313 -0.625,-26.25 -46.031,-0.075 -0.219,-38.363 -59.375,-0.937 -1.562,16.562 -33.125,-0.312 1.25,-27.5 -20.625,0 0.519,-21.489 -27.398,0 2.391,-45.605 17.937,0.527 0.883,-45.961 -49.996,0.231 -0.383,43.965 18.559,0 -1.766,47.726 -26.516,-0.883 -0.886,71.594 -53.914,-1.547 -0.575,-40.156 -55.1129,-1.605 -0.1563,20.636 -38.8867,-1.609 -0.5191,22.148 -10.34809,-0.625 z"
mpl_path = parse_path(bordersvg)
coords = mpl_path.to_polygons()
coords=rdp(coords[0], epsilon=3) #235 nodes 


svgpath="m 395.207,462.613 -7.734,0 -0.219,26.625 -161.531,-1.438 0.386,47.707 4.688,0.625 -0.313,-5.937 45,1.875 0.313,41.722 49.844,-0.629 0,5.782 13.281,0.312 0,-24.687 45.937,0 0.282,27.074 33.156,0.113 0.348,10.184 15.933,-0.809 0.906,-37.812 14.375,0 0,-7.188 -8.125,-0.625 1.563,-66.875 -47.836,-0.183 -0.254,-15.836 z"
mpl_path = parse_path(svgpath)
hole32 = mpl_path.to_polygons()
hole32=rdp(hole32[0], epsilon=3) 
#print(len(hole32)) #25
hole32=[tuple(i) for i in hole32]
holes.append([hole32,"32"])

svgpath="m 325.797,652.538 0.156,-9.687 -45,0 -0.156,17.812 -8.125,0.157 -0.156,-23.125 62.656,0.625 0.312,14.531 -9.687,-0.313 z"
mpl_path = parse_path(svgpath)
hole = mpl_path.to_polygons()
hole=rdp(hole[0], epsilon=3) 
#print(len(hole)) #9
hole=[tuple(i) for i in hole]
holes.append([hole,"16"])

svgpath="m 351.785,622.702 -0.312,-4.918 -26.371,-0.476 0.656,-32.965 20.172,0.133 0,-22.981 17.707,0.309 0.441,22.078 -4.836,0.024 -0.051,39.019 -7.406,-0.223 z"
mpl_path = parse_path(svgpath)
hole = mpl_path.to_polygons()
hole=rdp(hole[0], epsilon=3) 
#print(len(hole)) #9
hole=[tuple(i) for i in hole]
holes.append([hole,"18"])

svgpath="m 394.078,678.32 0,-22.5 55.938,0.156 -0.157,7.266 -29.843,-0.235 0,7.344 30,0.625 -0.157,8.906 -55.781,-1.562 z"
mpl_path = parse_path(svgpath)
hole = mpl_path.to_polygons()
hole=rdp(hole[0], epsilon=3) 
#print(len(hole)) #9
hole=[tuple(i) for i in hole]
holes.append([hole,"20"])

p="m 20.1559,421.296 0.2851,-9.914 12.493,0.32 0.1012,-4.222 17.6128,0.359 0.0668,14.004 -30.5589,-0.547 z"
mpl_path = parse_path(p)
hole = mpl_path.to_polygons()
hole=rdp(hole[0], epsilon=3) 
#print(len(hole)) #5
hole=[tuple(i) for i in hole]
holes.append([hole,"22"])

p="m 281.852,409.468 -44.196,-0.664 0.571,-17.328 -10.075,-8.746 -0.882,41.543 53.918,0.883 0.664,-15.688 z"
mpl_path = parse_path(p)
hole = mpl_path.to_polygons()
hole=rdp(hole[0], epsilon=3) 
#print(len(hole)) #7
hole=[tuple(i) for i in hole]
holes.append([hole,"24"])

p="m 56.902,458.081 -0.4411,-7.957 68.4731,0.496 -0.102,8.133 -67.93,-0.672 z"
mpl_path = parse_path(p)
hole = mpl_path.to_polygons()
hole=rdp(hole[0], epsilon=3) 
#print(len(hole)) #5
hole=[tuple(i) for i in hole]
holes.append([hole,"26"])

p="m 135.953,458.788 0,-7.187 26.094,0.469 0.312,6.562 -26.406,0.156 z"
mpl_path = parse_path(p)
hole = mpl_path.to_polygons()
hole=rdp(hole[0], epsilon=3) 
#print(len(hole)) #5
hole=[tuple(i) for i in hole]
holes.append([hole,"28"])

p="m 225.281,459.187 0.11,-7.184 65.531,0.739 -0.285,7.015 -65.356,-0.57 z"
mpl_path = parse_path(p)
hole = mpl_path.to_polygons()
hole=rdp(hole[0], epsilon=3) 
#print(len(hole)) #5
hole=[tuple(i) for i in hole]
holes.append([hole,"30"])

inners=[list(i[0]) for i in holes]
#print(type(inners[2][0]))
#print(inners[2])

border = Polygon(coords,holes=inners)
#border = Polygon(coords,holes=[[(1,10),(10,10),(10,1)],[(20,20),(30,30),(30,20)]])
#print(list(border.interiors[:]))

def graph():
    #plt.figure(figsize=(6, 12), dpi=100)
    '''
    xBorder, yBorder = zip(*coords)
    plt.gca().plot(xBorder, yBorder)
    '''
    for i in holes:
        xBorder, yBorder = zip(*i[0])
        plt.gca().plot(xBorder, yBorder, color="blue",linestyle=(0, (1, 1)))
        #plt.text(Polygon(i[0]).representative_point().coords[0][0], Polygon(i[0]).representative_point().coords[0][1], i[1])

    xBorder, yBorder = zip(*border.exterior.coords)
    plt.gca().plot(xBorder, yBorder,color="blue")
'''
    for i in border.buffer(-10):
        coordX = rdp(i.exterior.coords[:], epsilon=3) 
        if (len(i.interiors)>0):
            print(1)
            for ring in i.interiors:
                Px, Py = zip(*ring.coords)
                plt.gca().plot(Px, Py,color="black")

        Px, Py = zip(*coordX)
        plt.gca().plot(Px, Py,color="black")
        #plt.text(i.representative_point().coords[0][0], i.representative_point().coords[0][1], '({:0.2f}, {:0.2f})'.format(i.representative_point().coords[0][0],i.representative_point().coords[0][1]))

    plt.show() #polygons and stuff have attribute _ndim
'''
def graphinnermap():
    for i in border.buffer(-10):
        store=[]
        reduce(poly=i, epsilon=3, buffer=-10, areamin=1, store=store)

        for c in store:
            xBottom, yBottom = zip(*c)
            plt.gca().plot(xBottom, yBottom, color="black", linestyle='--', linewidth=1.0)  # marker=''

    coord1 = copy.deepcopy(border.exterior.coords[:])
    xBottom, yBottom = zip(*coord1)
    plt.gca().plot(xBottom, yBottom, color="black", linewidth=1.0)

    # plt.title("after cleaning")
    #plt.savefig('C:/Users/space/Documents/comp400 files/contour.png')
    plt.show()# almost look like a line [[310.746107311942, 430.23952726844936], [311.9021037673832, 415.66495666723335], [310.746107311942, 430.23952726844936]]

graph()
graphinnermap()