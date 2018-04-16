#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 21 16:56:53 2018

@author: duynguyen
Codes based on /Research/HiCDA/Packages/ttk/python/create-vtk-sparse.py
"""

import sys
sys.path = ['..']+sys.path

from pyvtk import *
import numpy as np
import math

def flat_list(list_of_list):
    "flat a list of list with unique elements"
    l = [item for sublist in list_of_list for item in sublist]
    return(list(set(l)))

def get_pixel_id(point_id, n, use_nbhd = False):
    """This creates a list of pixels containing the given
    point_id. The pixels' points are also returned"""
    xc = [ ]
    yc = [ ]
    if (use_nbhd):
        pixel_delta_x = [0, 1, 0, 1]
        pixel_delta_y = [0, 0, 1, 1]
    else :
        pixel_delta_x = [0]
        pixel_delta_y = [0]
    xp_ = point_id % n
    yp_ = math.floor(point_id / n)
    ## get pixels
    xc_ = xp_ - 1
    yc_ = yp_ - 1
    for i in range(len(pixel_delta_x)):
        xc_temp = xc_ + pixel_delta_x[i]
        yc_temp = yc_ + pixel_delta_y[i]
        if (xc_temp >= 0 and yc_temp >= 0 and xc_temp < (n-1) and yc_temp < (n-1)):
            xc.append(xc_temp)
            yc.append(yc_temp)
    pixels_select_old = [ xc[i] + yc[i]*(n-1) for i in range(len(xc))]
    return(pixels_select_old)

def get_point_id(pixel_id, n):
    """ This creates a list of the given pixel's 4 point ids.
    get_point_id(3,4) = [5,9,4,8]"""
    xc = pixel_id % (n-1)
    yc = math.floor(pixel_id / (n-1))
    xp = [xc+1, xc+1, xc, xc]
    yp = [yc, yc+1, yc, yc+1]
    return([n*y+x for x,y in zip(xp, yp)])

if __name__ == '__main__':
    ### * Inputs
    inputs = sys.argv[1:]
    filename = inputs[0] + "sim.vtk"
    filename_points = inputs[0] + "vtk-coords.csv"
    title = 'Unstructured Grid Example'

    #### * create height function
    filename_height = inputs[0] + 'heights-scalar.csv'
    heights = np.genfromtxt(filename_height, delimiter=',')
    ## delete colnames
    heights = np.delete(heights, (0), axis=0)

    ## dimension of hic matrix
    n = int(inputs[1])
    print('Dimension of HiC matrix: ', n)

    ## get default value for excluded d_height
    d_height_excluded = float(inputs[2])

    #### * sparse setting
    ## points_id_select: selection of points with id in square-coordinate (old coord.)
#    points_id_select = heights[:, 3] - 1
    points_id_select = heights[:, 1] - 1
    points_id_select = points_id_select.astype(int)    
    
    ## if python 2.7.x: use mydict.iteritems()
    pixels_select = [ get_pixel_id(point_id, n) for point_id in points_id_select]
    pixels_select = flat_list(pixels_select)

    pixels = [ get_point_id(pixel_id, n) for pixel_id in pixels_select]
    points_id_full_select = flat_list(pixels)

    ## create height function
#    d_height = heights[:, 2]
    d_height = heights[:, 0]
    heightDic = {p_id:d for p_id, d in zip(points_id_select, d_height) }
    points_id_diff = list(set(points_id_full_select) - set(points_id_select))
    heightDic.update({p_id:d_height_excluded for p_id in points_id_diff }  )

    ## convert points_id_full_select to 3D points
    px = np.array([ p_id % n for p_id in points_id_full_select])
    py = np.array([math.floor(p_id / n) for p_id in points_id_full_select])
    pz = np.repeat(0, len(px))
    points = np.column_stack((px, py, pz))
    
    ## output d_height for vtk
    d_height_vtk = [heightDic[p_id] for p_id in points_id_full_select]
    d_height_vtk = np.array(d_height_vtk)
    
    ## convert pixels to pixels in vtk coord (i.e. in order of points)
    map_point_id_square_to_vtk_coord = {old:new for new, old in enumerate(points_id_full_select)}
    pixels_vtk_coord = [ ]
    for row in pixels:
        pixels_vtk_coord.append([ map_point_id_square_to_vtk_coord[val] for val in row])

    #### * write to files
    vtk = VtkData(UnstructuredGrid(points = points, pixel = pixels_vtk_coord),
                  PointData(Scalars( d_height_vtk, name = 'heights' )),
                  'Unstructured Grid Example')    
    np.savetxt(filename_points, points[:,[0,1]], delimiter=",")    
    vtk.tofile(filename)
    
    #vtk.tofile(filename, 'binary') # cause some errors
