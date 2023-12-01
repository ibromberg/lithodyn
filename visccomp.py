#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 24 16:40:00 2023

@author: ibromberg

viscosity diffs between 17, 16.5, and 16
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from matplotlib import cm

def cart2sph(x,y,z): # function for converting cartesian to spherical
    theta = np.arctan2(y,x) # azimuth
    phi = np.arctan2(z,np.sqrt(x**2 + y**2)) # elevation
    r = np.sqrt(x**2 + y**2 + z**2)#-6371000
    return theta, phi, r

def sph2cart(azimuth,elevation,r): # function for converting spherical to cartesian
    x = r * np.cos(elevation) * np.cos(azimuth)
    y = r * np.cos(elevation) * np.sin(azimuth)
    z = r * np.sin(elevation)
    return x, y, z

def process(df): # function to pull x, y, and z out of a dataframe and convert it to r theta phi
    x = np.array(df["lat"].values.tolist())
    y = np.array(df["long"].values.tolist())
    z = np.array(df["visc"].values.tolist())
    theta, phi, r = cart2sph(x,y,z)
    
    rsub=[]
    for i in range(0,len(r),1):
        rsub.append(r[i]-6371000)
    
    return x,y,z

############

t16 = pd.read_csv('viscosity_log10__16.0.out',
                  sep="       ",names=["lat","long","visc"])
t165 = pd.read_csv('viscosity_log10__16.5.out',
                   sep="       ",names=["lat","long","visc"])
t17 = pd.read_csv('viscosity_log10__17.0.out',
                  sep="       ",names=["lat","long","visc"])


# read in the states data
arizona = pd.read_csv('deformed_state_boundaries/arizona_17.0.dat', sep = ' ', header=None)
california = pd.read_csv('deformed_state_boundaries/california_17.0.dat', sep = ' ', header=None)
colorado = pd.read_csv('deformed_state_boundaries/colorado_17.0.dat', sep = ' ', header=None)
mexico = pd.read_csv('deformed_state_boundaries/mexico_17.0.dat', sep = ' ', header=None)
nevada = pd.read_csv('deformed_state_boundaries/nevada_17.0.dat', sep = ' ', header=None)
new_mexico = pd.read_csv('deformed_state_boundaries/new_mexico_17.0.dat', sep = ' ', header=None)
oregon = pd.read_csv('deformed_state_boundaries/oregon_17.0.dat', sep = ' ', header=None)
texas = pd.read_csv('deformed_state_boundaries/texas_17.0.dat', sep = ' ', header=None)
utah = pd.read_csv('deformed_state_boundaries/utah_17.0.dat', sep = ' ', header=None)
wyoming = pd.read_csv('deformed_state_boundaries/wyoming_17.0.dat', sep = ' ', header=None)

theta16, phi16, v16 = (np.array(t16["lat"].values.tolist()), 
                       np.array(t16["long"].values.tolist()), 
                       np.array(t16["visc"].values.tolist()))
theta165, phi165, v165 = (np.array(t165["lat"].values.tolist()), 
                          np.array(t165["long"].values.tolist()), 
                          np.array(t165["visc"].values.tolist()))
theta17, phi17, v17 = (np.array(t17["lat"].values.tolist()), 
                       np.array(t17["long"].values.tolist()), 
                       np.array(t17["visc"].values.tolist()))

f, ax = plt.subplots(3)

f.set_size_inches(10, 15)
f.suptitle("Comparison of visc")

ax[0].set_aspect('equal')
ax[1].set_aspect('equal')
ax[2].set_aspect('equal')

m=ax[0].tricontourf(theta16,phi16,v16)
ax[0].set_title("16Ma")

ax[1].tricontourf(theta165,phi165,v165)
ax[1].set_title("16.5Ma")

ax[2].tricontourf(theta17,phi17,v17)
ax[2].set_title("17Ma")

# one color bar for them to share
f.subplots_adjust(right=0.8)
cbar_ax = f.add_axes([0.85,0.15,0.03,0.7])
f.colorbar(m, cax=cbar_ax)

#plt.savefig("viscositycomp.png")

# just 17 Ma
f, ax = plt.subplots(1)
f.set_size_inches(10, 5)

plt.tricontourf(theta17,phi17,v17,levels=100)
plt.xlim([np.min(theta17),np.max(theta17)])
plt.ylim([np.min(phi17),np.max(phi17)])
plt.colorbar(label="$log_{10}$")
plt.title("Vertically Averaged Effective Viscosity at 17Ma")

# add the deformed state lines
plt.plot(arizona[0],arizona[1], color = "black")
plt.plot(california[0],california[1], color = "black")
plt.plot(colorado[0],colorado[1], color = "black")
#plt.plot(mexico[0],mexico[1], color = "black")
plt.plot(nevada[0],nevada[1], color = "black")
plt.plot(new_mexico[0],new_mexico[1], color = "black")
plt.plot(oregon[0],oregon[1], color = "black")
plt.plot(texas[0],texas[1], color = "black")
plt.plot(utah[0],utah[1], color = "black")
plt.plot(wyoming[0],wyoming[1], color = "black")

plt.savefig("poster/visco.png",dpi=300,bbox_inches='tight')