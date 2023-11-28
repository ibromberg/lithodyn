#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 12:19:51 2023

@author: ibromberg

creating my own differential velocity magnitude file

theta is long, phi is lat
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import cm

def cart2sph(x,y,z): # function for converting cartesian coords to spherical
    theta = np.arctan2(y,x) * 180/np.pi # azimuth
    phi = np.arctan2(z,np.sqrt(x**2 + y**2)) * 180/np.pi # elevation
    r = np.sqrt(x**2 + y**2 + z**2)#-6371000
    return theta, phi, r


def cart2sphV(long,lat,vx,vy,vz): # convert velocities to spherical
    theta = long * np.pi/180
    phi = [90-(i * np.pi/180) for i in lat] 

    #vr = vx * np.sin(phi) * np.cos(theta) + vy * np.sin(phi) * np.sin(theta) + vz * np.cos(phi)
    vphi = vx * np.cos(phi) * np.cos(theta) + vy * np.cos(phi) * np.sin(theta) - vz * np.sin(phi)
    vtheta = -vx * np.sin(theta) + vy * np.cos(theta)
    
    vLong = vtheta
    vLat = -vphi
    
    return vLong, vLat

def diff(vTopo,v28):
    vDiff=[]
    for i in range(len(vTopo)):
        vDiff.append(v28[i]-vTopo[i])
    return vDiff

# -------------------------------------------------------

# read in data from comsol
topo = pd.read_csv('fcm5_lc22_isostacy_noslip_bottom.txt',sep=" ",header=None,names=["x","y","z","elev"])
vtopo = pd.read_csv('fcm2_nobc_isotest_wallnoslip_bottomconstraint_topov.txt',sep=" ",comment='%',header=None,names=["x","y","z","vx","vy","vz"])
#vlc = pd.read_csv('fcm5_lc203_iso_noslip_bottom_lcv.txt',sep=" ",comment='%',header=None,names=["x","y","z","vx","vy","vz"])
v28 = pd.read_csv('fcm2_nobc_isotest_wallnoslip_bottomconstraint_28kmv.txt',sep=" ",comment='%',header=None,names=["x","y","z","vx","vy","vz"])
#vlc = pd.read_csv('fcm5_lc22_isostacy_noslip_bottom.txt',sep=" ",header=None,names=["x","y","z","vx","vy","vz"])

#v28 = v28.iloc[::20, :]

# arrays for surface topography
x, y, z = np.array(topo["x"].values.tolist()), np.array(topo["y"].values.tolist()), np.array(topo["z"].values.tolist())
theta, phi, r = cart2sph(x,y,z)

# arrays for surface XYZ coords and velocities in XYZ 
xsurf, ysurf, zsurf = np.array(vtopo["x"].values.tolist()), np.array(vtopo["y"].values.tolist()), np.array(vtopo["z"].values.tolist())
thetasurf, phisurf, rsurf = cart2sph(xsurf,ysurf,zsurf) # convert surface to theta and phi


vxsurf, vysurf, vzsurf = np.array(vtopo["vx"].values.tolist()), np.array(vtopo["vy"].values.tolist()), np.array(vtopo["vz"].values.tolist())
vLongsurf, vLatsurf = cart2sphV(thetasurf,phisurf,vxsurf,vysurf,vzsurf)
#vLongsurf = np.array([1 for i in vLongsurf])
#vLatsurf = np.array([1 for i in vLatsurf])

# arrays for LC xyz coords and velocities in xyz
"""
xlc, ylc, zlc = np.array(vlc["x"].values.tolist()), np.array(vlc["y"].values.tolist()), np.array(vlc["z"].values.tolist())
thetalc, philc, rlc = cart2sph(xlc,ylc,zlc) # convert surface to theta and phi

vxlc, vylc, vzlc = np.array(vlc["vx"].values.tolist()), np.array(vlc["vy"].values.tolist()), np.array(vlc["vz"].values.tolist())
vLonglc, vLatlc = cart2sphV(thetalc,philc,vxlc,vylc,vzlc)
"""

# arrays for 28 km below xyz coords and velocities in xyz
x28, y28, z28 = np.array(v28["x"].values.tolist()), np.array(v28["y"].values.tolist()), np.array(v28["z"].values.tolist())
theta28, phi28, r28 = cart2sph(x28,y28,z28) # convert surface to theta and phi

vx28, vy28, vz28 = np.array(v28["vx"].values.tolist()), np.array(v28["vy"].values.tolist()), np.array(v28["vz"].values.tolist())
vLong28, vLat28 = cart2sphV(theta28,phi28,vx28,vy28,vz28)


# calculate differential flow

vLongdiff = diff(vLongsurf,vLong28)
vLatdiff = diff(vLatsurf,vLat28)

# -------------------------------------------------------

f, ax = plt.subplots(3)
f.set_size_inches(10, 15)
f.suptitle("Differential Flow: FCM5, LC 20.3")

ax[0].set_aspect('equal')
ax[1].set_aspect('equal')
ax[2].set_aspect('equal')

im0=ax[0].tricontourf(theta,phi,r,cmap="gist_earth") # topography
ax[0].quiver(thetasurf,phisurf, vLongsurf, vLatsurf) # surface arrows
ax[0].set_title("Surface Flow")

"""
im1=ax[1].tricontourf(theta,phi,r,cmap="gist_earth") # topography
ax[1].quiver(thetalc,philc, vLonglc, vLatlc) # arrows
ax[1].set_title("LC Flow")
"""

im1=ax[1].tricontourf(theta,phi,r,cmap="gist_earth") # topography
ax[1].quiver(theta28,phi28, vLong28, vLat28) # surface arrows
ax[1].set_title("28 Km Below Surface")

im2=ax[2].tricontourf(theta,phi,r,cmap="gist_earth") # topography
ax[2].quiver(thetasurf,phisurf, vLongdiff, vLatdiff) # arrows
ax[2].set_title("Differential Flow")

plt.savefig("diffflow.png",bbox_inches='tight')
