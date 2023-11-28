#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Nov 13 16:21:41 2023

@author: ibromberg

create grid of points for comsol to export velocities at
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# theta = long, phi = lat

def cart2sph(x,y,z): # function for converting cartesian to spherical
    azimuth = np.arctan2(y,x) * 180/np.pi #theta, long
    elevation = np.arctan2(z,np.sqrt(x**2 + y**2)) * 180/np.pi #phi, lat
    r = np.sqrt(x**2 + y**2 + z**2)
    return elevation, azimuth, r

def sph2cart(long,lat,r): # function for converting spherical to cartesian
    long = long*np.pi/180
    lat = lat*np.pi/180
    x = r * np.cos(long) * np.cos(lat)
    y = r * np.cos(long) * np.sin(lat)
    z = r * np.sin(long)
    return x, y, z

topo = pd.read_csv('fcm2_nobc_isotest_wallnoslip_bottomconstraint_topo.txt',
                   sep=" ", comment='%',
                   header=None,names=["x","y","z","elev"])

# make data set to work with smaller, every 10th value
#topo = topo.iloc[::10, :]

# write x, y, and z back into a file for topography
x, y, z = (np.array(topo["x"].values.tolist()), 
                np.array(topo["y"].values.tolist()), 
                np.array(topo["z"].values.tolist()))
np.savetxt('datagrid_topo.txt',np.transpose([x,y,z]))

# make lat long r for 28km bsl
lat, long, r = cart2sph(x,y,z)
r = np.array([6343000 for i in r]) # make all r uniform

# write x,y,z into file for 28km bsl
x, y, z = sph2cart(lat, long, r)
np.savetxt('datagrid_28kmBSL.txt', np.transpose([x,y,z]))