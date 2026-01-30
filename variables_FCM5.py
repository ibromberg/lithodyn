#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Dec  1 15:28:11 2025

@author: ibromberg
"""

import numpy as np
from scipy.linalg import eig
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import pandas as pd
import scipy.optimize

# fcm 5

lc =            [10,     15,     20,     25,     30,     50,     90]

diff =          [934.28, 823.24, 753.75, 705.03, 668.04, 577.55, 492.29] #ali max - 16Ma max
# these are all at 17Ma
normalsurfv =   [5.1,    4.44,   4.05,   3.77,   3.56,   3.03,   2.56] # u*nx + v*ny + w*nz
surfv =         [5.73,   5.08,   4.69,   4.42,   4.21,   3.72,   3.28]   # spf.U
slicev =        [8.4,    6.8,    5.97,   5.44,   5.07,   4.24,   3.51] # using 3D slices (in LC)
tanv =          [6.19,   4.79,   4.07,   3.90,    3.76,   3.39,   2.96]
minvisc =       [20.3,   20.5,   20.6,   20.7,   20.8,   21,     21] # from comsol

# tanv = []
# for i in range(0,len(diff)):
#     vtan = np.sqrt(surfv[i]**2 - normalsurfv[i]**2)
#     tanv.append(vtan)

#lcvisc = [10, 15, 20, 25, 30, 40, 50, 60, 70, 80, 90] # from matlab
#minvisc = [21.23, ]


lc2 = [10, 30, 50]

fcm3_diff = [990.97, 693.36, 594.83]
fcm3_normsurfv = [5.7, 3.95, 3.32]
fcm3_tanv = [6.53, 4.24, 3.74]
fcm3_slicev = [9.02, 5.48, 4.57]
fcm3_surfv = [6.6, 4.8, 4.19]

fcm7_diff = [910.07, 659.08, 572.78]
fcm7_normsurfv = [4.84, 3.4, 2.89]
fcm7_tanv = [6.04, 3.52, 3.21]
fcm7_slicev = [8.16, 4.9, 4.1]
fcm7_surfv = [5.34, 3.95, 3.46]

savepath = "/Users/ibromberg/Documents/COMSOL 61/extrusion/plots/variables.pdf"
savepath = "/Users/ibromberg/Documents/COMSOL 63/prelim/variables.pdf"

#-------------

fig, ax = plt.subplots(2, 3, dpi=300)
fig.set_size_inches(13, 10)

ax[0,0].plot(lc2,fcm3_diff,marker='o')
ax[0,0].plot(lc,diff,marker='*')
ax[0,0].plot(lc2,fcm7_diff,marker='s')
ax[0,0].set_xlabel("Percentage LC Strength")
ax[0,0].set_ylabel("Elevation Difference (m)")
ax[0,0].set_title("Difference between Bahadori and simulation \n at central elevation peak at 16Ma")

ax[0,1].plot(lc2,fcm3_normsurfv,marker='o')
ax[0,1].plot(lc,normalsurfv,marker='*')
ax[0,1].plot(lc2,fcm7_normsurfv,marker='s')
ax[0,1].set_xlabel("Percentage LC Strength")
ax[0,1].set_ylabel("Velocity (mm/yr)")
ax[0,1].set_title("Max Normal Surface Velocity at 17Ma")
ax[0,1].legend(labels=["FCM3", "FCM5", "FCM7"])

ax[0,2].plot(lc2,fcm3_surfv,marker='o')
ax[0,2].plot(lc,surfv,marker='*')
ax[0,2].plot(lc2,fcm7_surfv,marker='s')
ax[0,2].set_xlabel("Percentage LC Strength")
ax[0,2].set_ylabel("Velocity (mm/yr)")
ax[0,2].set_title("Max Surface Velocity Magnitude at 17Ma")

ax[1,0].plot(lc2,fcm3_slicev,marker='o')
ax[1,0].plot(lc,slicev,marker='*')
ax[1,0].plot(lc2,fcm7_slicev,marker='s')
ax[1,0].set_xlabel("Percentage LC Strength")
ax[1,0].set_ylabel("Velocity (mm/yr)")
ax[1,0].set_title("Max Velocity in LC at 17Ma")

ax[1,1].plot(lc2,fcm3_tanv,marker='o')
ax[1,1].plot(lc,tanv,marker='*')
ax[1,1].plot(lc2,fcm7_tanv,marker='s')
ax[1,1].set_xlabel("Percentage LC Strength")
ax[1,1].set_ylabel("Velocity (mm/yr)")
ax[1,1].set_title("Max Tangential Surface Velocity at 17Ma")

ax[1,2].plot(lc,minvisc,color='k',marker='.')
ax[1,2].set_xlabel("Percentage LC Strength")
ax[1,2].set_ylabel("Viscosity Magnitude")
ax[1,2].set_title("Minimum LC Viscosity")

plt.show()

fig, ax = plt.subplots(1,2, dpi=300)
fig.set_size_inches(10, 5)


ax[0].plot(lc2,fcm3_normsurfv,marker='o')
ax[0].scatter(50,3.26,color='k')
ax[0].plot(lc,normalsurfv,marker='*')
ax[0].plot(lc2,fcm7_normsurfv,marker='s')
ax[0].set_xlabel("Percentage LC Strength")
ax[0].set_ylabel("Velocity (mm/yr)")
ax[0].set_title("Max Normal Surface Velocity at 17Ma")
ax[0].legend(labels=["FCM3", "FCM4", "FCM5", "FCM7"])

ax[1].plot(lc2,fcm3_tanv,marker='o')
ax[1].plot(lc,tanv,marker='*')
ax[1].plot(lc2,fcm7_tanv,marker='s')
ax[1].set_xlabel("Percentage LC Strength")
ax[1].set_ylabel("Velocity (mm/yr)")
ax[1].set_title("Max $\sqrt{ v_{lat}^2 + v_{long}^2 }$ Surface Velocity at 17Ma")

#plt.savefig(savepath,bbox_inches='tight')
plt.show()

#######################
lc = [10, 20, 30, 50, 90]
fcm3 = [5.7, 4.5, 3.95, 3.32, 2.75]
fcm5 =   [5.1, 4.05, 3.56, 3.03, 2.56] # u*nx + v*ny + w*nz
fcm7 = [4.84, 3.85, 3.4, 2.89, 2.46] # (up)

fig, ax = plt.subplots()
fig.set_size_inches(5,5)
plt.scatter(lc, fcm3, marker='o',label="FCM 3")
plt.scatter(lc, fcm5, marker='*',label="FCM 5")
plt.scatter(lc, fcm7, marker='s',label="FCM 7")
plt.title("Max Normal Surface Velocity")
plt.ylabel("Velocity (mm/yr)")
plt.xlabel("Percentage Lower Crust Strength")
plt.legend()

# for lc 21
# fcm 5: 3.68

xs = np.linspace(0,100,100)

def monoExp(x,m,t,b):
    return m*np.exp(-t*x)+b

def fit(fcm):
    params, cv = scipy.optimize.curve_fit(monoExp, lc, fcm)
    m,t,b = params
    return monoExp(xs,m,t,b)


plt.plot(xs,fit(fcm3))
plt.plot(xs,fit(fcm5))
plt.plot(xs,fit(fcm7))
plt.savefig(savepath,bbox_inches='tight')
