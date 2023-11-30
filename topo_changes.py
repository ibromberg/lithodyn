#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Nov 30 13:51:00 2023

@author: ibromberg

just plot some topography elevation changes
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

def cart2sph(x,y,z): # function for converting cartesian coords to spherical
    theta = np.arctan2(y,x) * 180/np.pi # azimuth
    phi = np.arctan2(z,np.sqrt(x**2 + y**2)) * 180/np.pi # elevation
    r = np.sqrt(x**2 + y**2 + z**2)-6371000
    return theta, phi, r

def process(df):
    x = np.array(df["x"].values.tolist())
    y = np.array(df["y"].values.tolist())
    z = np.array(df["z"].values.tolist())
    
    theta, phi, r = cart2sph(x,y,z)
    
    return theta, phi, r

############
# open files

start = pd.read_csv('fcm2_nobc_isotest_wallnoslip_bottomconstraint_topo_0.txt',
                    sep=" ", comment='%', header=None,
                    names=["x","y","z","elev"])
mid = pd.read_csv('fcm2_nobc_isotest_wallnoslip_bottomconstraint_topo_5e5.txt',
                    sep=" ", comment='%', header=None,
                    names=["x","y","z","elev"])
final = pd.read_csv('fcm2_nobc_isotest_wallnoslip_bottomconstraint_topo_1e6.txt',
                    sep=" ", comment='%', header=None,
                    names=["x","y","z","elev"])


########

theta_0, phi_0, r_0 = process(start)
theta_5e5, phi_5e5, r_5e5 = process(mid)
theta_1e6, phi_1e6, r_1e6 = process(final)


############ PLOT
"""
f, ax = plt.subplots(1)
#f.set_size_inches(12, 10)
levels = np.linspace(-5000,5000,51)

# surface and 28km bsl velocities
im2=plt.tricontourf(theta_0,phi_0,r_0,levels=levels,cmap="gist_earth") # topography

#plt.xlim([-122.5,-105])
#plt.ylim([28,42])
plt.title("Elevation at t=0")
f.colorbar(im2)"""

fig, ax = plt.subplots(3)
fig.set_size_inches(5, 10)
fig.suptitle("Elevation Through Time")
fig.tight_layout()

ax[0].set_aspect('equal')
ax[1].set_aspect('equal')
ax[2].set_aspect('equal')

levels = 100#np.linspace(-5000,5000,51)

im0=ax[0].tricontourf(theta_0,phi_0,r_0,levels=levels,cmap="gist_earth") # without bottom constraint
ax[0].set_title("t = 0")
#fig.colorbar(im0,ax=ax[0])

im1=ax[1].tricontourf(theta_5e5,phi_5e5,r_5e5,levels=levels,cmap="gist_earth") # with bottom constraint
ax[1].set_title("t = 5e5")
#fig.colorbar(im1,ax=ax[1])

im2=ax[2].tricontourf(theta_1e6,phi_1e6,r_1e6,levels=levels,cmap="gist_earth") # difference between the two
ax[2].set_title("t = 1e6")
#fig.colorbar(im2,ax=ax[2])

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85,0.15,0.03,0.7])
fig.colorbar(im0, cax=cbar_ax)

plt.savefig("poster/elevation.png",dpi=300,bbox_inches='tight')