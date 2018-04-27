#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jan 29 13:26:16 2018

@author: duynguyen
Codes based on /Research/HiCDA/Packages/ttk/python/write-vtk-upper-triangle.py
"""

import sys
sys.path = ['..']+sys.path

sys.path.append("/p/keles/DBChIP/volume4/HiCDA/Packages/pyvtk")

from pyvtk import *
import numpy as np

if __name__ == '__main__':
    ### Inputs
    inputs = sys.argv[1:]
    ## tested on n = 5000
    #filename = '/home/duynguyen/ResearchLocal/HiCDA/test-TreeHiC/temp/sim.vtk'
    filename = inputs[0] + "sim.vtk"
    title = 'Unstructured Grid Example'
    pixels = []
    points = []
    triangles = []

    #### create height function
    filename_height = inputs[0] + 'heights-scalar.csv'
    #heights = np.genfromtxt('/home/duynguyen/ResearchLocal/HiCDA/test-TreeHiC/temp/heights-scalar.csv', delimiter=',')
    heights = np.genfromtxt(filename_height, delimiter=',')
    heights = np.delete(heights, [0])

    n = int((-3 + (17 + 8*len(heights)) ** 0.5 )/2)
    print('Dimension of HiC matrix: ', n)
    #### create point data
    px = np.array(list(range(0,n)))
    px = np.tile(px, n)
    py = np.array(list(range(0, n)))
    py = np.repeat(py, n, axis=0)
    pz = np.repeat(0,n*n)
    points = np.column_stack((px, py, pz))

    # @
    # upTriId = [0,4,8,12,1,5,9,13,6,10,14,11,15] # n=4
    upTriId = np.arange(0, n*n, step = n)
    for i in range(1, n):
        temp = upTriId + 1
        upTriId = np.append(upTriId, temp[-(n- i + 1):])

    points = points[upTriId]
    mydict = {old:new for new, old in enumerate(upTriId)}

    #### create cell data
    temp = np.array([1, n+1,0 ,n])
    #temp = [1, n+1,0 ,n]
    for block in range(0,n-1):
        for j in range(1,n):
            pixels.append(temp + j - 1)
        temp = pixels[len(pixels) - 1] + 2

    pixels = np.array(pixels)
    # upTriCellId = [0,3,6,4,7,8] # n=4
    upTriCellId = np.arange(0, (n-1)*(n-1), step = n-1)
    for i in range(1, n-1):
        temp = upTriCellId + 1
        upTriCellId = np.append(upTriCellId, temp[-(n- i -1):])


    pixels = pixels[upTriCellId]
    pixels_new = [ ]
    for row in pixels:
        pixels_new.append([mydict[val] for val in row])


    #### write to file
    vtk = VtkData(UnstructuredGrid(points = points, pixel = pixels_new),
            PointData(Scalars( heights, name = 'heights' )),
            'Unstructured Grid Example')

    vtk.tofile(filename)
    #vtk.tofile(filename, 'binary') # cause some errors
    #print(open(filename).read())
