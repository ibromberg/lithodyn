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
import matplotlib.patheffects as pe

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
        vDiff.append(vTopo[i]-v28[i])
    return vDiff

def draw_line(lat,long,angle,length):
  cartesianAngleRadians = (450-angle)*np.pi/180.0
  long_end = long + length * np.cos(cartesianAngleRadians)
  lat_end = lat + length * np.sin(cartesianAngleRadians)
  #print(str(long) + " " + str(long_end))
  #print(str(lat) + " " + str(lat_end))
  #plt.plot([long, long_end],[lat,lat_end],color = "magenta",linewidth=3)
  return long_end, lat_end
  

# -------------------------------------------------------

# read in data from comsol: topography, velocity topo, and velocity 28kmbsl
topo = pd.read_csv('fcm5_lc22_isostacy_noslip_bottom.txt',sep=" ",
                   header=None,names=["x","y","z","elev"])
vtopo = pd.read_csv('fcm2_nobc_isotest_wallnoslip_bottomconstraint_topov.txt',
                    sep=" ",comment='%',header=None,names=["x","y","z","vx","vy","vz"])
v28 = pd.read_csv('fcm2_nobc_isotest_wallnoslip_bottomconstraint_28kmv.txt',
                  sep=" ",comment='%',header=None,names=["x","y","z","vx","vy","vz"])

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

# read in core complex data
core = pd.read_csv('corecomplex.txt', sep=" ", header=None, comment='%', 
                   names=['name','long','lat','deg'])

# ------------------------------------------

# arrays for surface topography
x, y, z = (np.array(topo["x"].values.tolist()), 
           np.array(topo["y"].values.tolist()), 
           np.array(topo["z"].values.tolist()))
theta, phi, r = cart2sph(x,y,z)

# arrays for surface XYZ coords
xsurf, ysurf, zsurf = (np.array(vtopo["x"].values.tolist()), 
                       np.array(vtopo["y"].values.tolist()), 
                       np.array(vtopo["z"].values.tolist()))
thetasurf, phisurf, rsurf = cart2sph(xsurf,ysurf,zsurf) # convert surface to theta and phi

# arrays for surface xyz velocities
vxsurf, vysurf, vzsurf = (np.array(vtopo["vx"].values.tolist()), 
                          np.array(vtopo["vy"].values.tolist()), 
                          np.array(vtopo["vz"].values.tolist()))
vLongsurf, vLatsurf = cart2sphV(thetasurf,phisurf,vxsurf,vysurf,vzsurf)

# arrays for 28 km below xyz coords and velocities in xyz
x28, y28, z28 = (np.array(v28["x"].values.tolist()), 
                 np.array(v28["y"].values.tolist()), 
                 np.array(v28["z"].values.tolist()))
theta28, phi28, r28 = cart2sph(x28,y28,z28) 

vx28, vy28, vz28 = (np.array(v28["vx"].values.tolist()), 
                    np.array(v28["vy"].values.tolist()), 
                    np.array(v28["vz"].values.tolist()))
vLong28, vLat28 = cart2sphV(theta28,phi28,vx28,vy28,vz28)


# calculate differential flow

vLongdiff = diff(vLongsurf,vLong28)
vLatdiff = diff(vLatsurf,vLat28)

# core complexes
cclong = core['long']#.values.tolist()
cclat = core['lat']#.values.tolist()
ccdeg = core['deg']

# -------------------------------------------------------

f, ax = plt.subplots(1)
f.set_size_inches(12, 10)

# surface and 28km bsl velocities
#plt.tricontourf(theta,phi,r,cmap="gist_earth") # topography
plt.quiver(thetasurf,phisurf, vLongsurf, vLatsurf,color='royalblue',label="Surface Velocity") # surface arrows
plt.quiver(theta28,phi28, vLong28, vLat28,color='firebrick',label="28 km BSL Velocity") # surface arrows

plt.xlim([-122.5,-105])
plt.ylim([28,42])
plt.title("Surface & 28km Below Sea Level Flow")

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

# plot core complexes
# dots = plt.scatter(mcc[0],mcc[1],marker='o',s=150,color="grey",edgecolors="black",linewidth=1,zorder=100)
plotcc = plt.scatter(cclong,cclat,marker='o',s=100,color='indigo')

for i in range(0,len(ccdeg)):
    #draw_line(cclong[i],cclat[i],ccdeg[i],10)
    long_end, lat_end = draw_line(cclat[i],cclong[i],ccdeg[i],1)
    plt.plot([cclong[i], long_end],[cclat[i],lat_end],color = "indigo",linewidth=3)

plt.savefig("poster/surfand28v.png",dpi=300,bbox_inches='tight')

# differential velocities
plt.clf()
plt.quiver(thetasurf,phisurf, vLongdiff, vLatdiff,color='royalblue',label="Differential Velocity") # surface arrows

plt.xlim([-122.5,-105])
plt.ylim([28,42])
plt.title("Differential Velocity")

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

# plot core complexes
# dots = plt.scatter(mcc[0],mcc[1],marker='o',s=150,color="grey",edgecolors="black",linewidth=1,zorder=100)
plotcc = plt.scatter(cclong,cclat,marker='o',s=100,color='indigo')

for i in range(0,len(ccdeg)):
    #draw_line(cclong[i],cclat[i],ccdeg[i],10)
    long_end, lat_end = draw_line(cclat[i],cclong[i],ccdeg[i],1)
    plt.plot([cclong[i], long_end],[cclat[i],lat_end],color = "indigo",linewidth=3)

plt.savefig("poster/diffv.png",dpi=300,bbox_inches='tight')