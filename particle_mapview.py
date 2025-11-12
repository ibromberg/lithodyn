#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Sep 16 12:01:37 2025

@author: ibromberg


particle tracing tracker
"""


import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from matplotlib import colors

from matplotlib import cm

# theta long phi lat

def cart2sph(x,y,z): # function for converting cartesian coords to spherical
    theta = np.arctan2(y,x) * 180/np.pi# azimuth
    phi = np.arctan2(z,np.sqrt(x**2 + y**2)) * 180/np.pi # elevation
    r = np.sqrt(x**2 + y**2 + z**2)-6371000
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

def norm(vLong,vLat):
    N = np.sqrt(vLong**2 + vLat**2)
    vLong = vLong/N 
    vLat = vLat/N 
    return vLong, vLat

def process(df):
    x = np.array(df["x"].values.tolist())
    y = np.array(df["y"].values.tolist())
    z = np.array(df["z"].values.tolist())
    
    theta, phi, r = cart2sph(x, y, z)
    
    return theta, phi, r

# -------------------------------------------------------

# imports

path1 = "/Users/ibromberg/Documents/COMSOL 63/particle data/fcm5lc203_17ma.txt"
path2 = "/Users/ibromberg/Documents/COMSOL 63/particle data/fcm5lc203_16ma.txt"

savepath = "/Users/ibromberg/Documents/COMSOL 63/plots/"

ma17 = pd.read_csv(path1, sep=" ",comment='%',header=None,
                    names=["index","x","y","z","vx","vy","vz"])
ma16 = pd.read_csv(path2, sep=" ",comment='%',header=None,
                    names=["index","x","y","z","vx","vy","vz"])

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

# core complexes
cclong = core['long']#.values.tolist()
cclat = core['lat']#.values.tolist()
ccdeg = core['deg']


# -------------------------------------------------------

#  processing

# southwest: -110 to -118, below 30 to 37
# core complexes: -111,-116, 33 to 36

rmin = -20000
rmax = 0
longmin = -116
longmax = -111
latmin = 33
latmax = 36

long_init, lat_init, r_init = process(ma17)
long_fin, lat_fin, r_fin = process(ma16)

ma16["lat"] = lat_init
ma16["long"] = long_init
ma16["r"] = r_init
ma17["lat"] = lat_fin
ma17["long"] = long_fin
ma17["r"] = r_fin


ma16 = ma16[(ma17.lat < latmax) & (ma17.long < longmax) & (ma17.long > longmin) & (ma17.r > rmin) & (ma17.r < rmax)]
ma17 = ma17[(ma17.lat < latmax) & (ma17.long < longmax) & (ma17.long > longmin) & (ma17.r > rmin) & (ma17.r < rmax)]

lat_init, long_init, r_init = (np.array( ma17["lat"].values.tolist() ),
                               np.array( ma17["long"].values.tolist() ),
                               np.array( ma17["r"].values.tolist() ))
lat_fin, long_fin, r_fin = (np.array( ma16["lat"].values.tolist() ),
                            np.array( ma16["long"].values.tolist() ),
                            np.array( ma16["r"].values.tolist() ))

r_diff = (r_fin - r_init)
lat_diff = (lat_fin-lat_init)
long_diff = (long_fin-long_init)

medr = np.median(r_init)

print(np.isnan(np.min(long_diff)))


# -------------------------------------------------------

# plotting
fig, ax = plt.subplots()

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

plt.xlim(longmin,longmax)
plt.ylim(latmin,latmax)

#plt.axvline(x = -112, color = 'black',linestyle='dashed')
#plt.axvline(x = -112.5, color = 'black',linestyle='dashed')

# set color bar 0 at 0
normalcolor = colors.TwoSlopeNorm(vmin=min(r_diff),vcenter=0,vmax=max(r_diff)) 


#plt.scatter(long_init,lat_init,c=r_diff, ec='black', linewidth=0.5, cmap='bwr', s=200, norm=normalcolor)
plt.tricontourf(long_init,lat_init,r_diff,cmap='bwr',norm=normalcolor)
#for n in range(0,len(r_init)):
#    plt.arrow(long_init[n],lat_init[n],long_diff[n],lat_diff[n])

scaling = 0.01 #5e-9 # change as needed; 0.01 
# 1 mm/yr = 3.171e-11 m/s
# 1 degree = 111,111 m; 1000 m = 111.111 deg
plt.subplots_adjust(bottom=0.2)

Q1 = plt.quiver(long_init,lat_init, long_diff, lat_diff,color='black',scale=scaling)
qk = ax.quiverkey(Q1, X=0.25, Y=-0.2, U=1/111111 * 1000, label='1000 m Horizontal Displacement',labelpos='E')

plt.colorbar(label="Vertical Particle Displacement (m)")

plt.title("Particle Displacement from 16Ma - 17Ma, Median Particle Start Depth " + str(round(abs(medr)/1000,1)) + "km")


# core complexes
plotcc = plt.scatter(cclong,cclat,marker='o',s=100,color='indigo',label="Core Complexes")

for i in range(0,len(ccdeg)):
    #draw_line(cclong[i],cclat[i],ccdeg[i],10)
    long_end, lat_end = draw_line(cclat[i],cclong[i],ccdeg[i],1)
    plt.plot([cclong[i], long_end],[cclat[i],lat_end],color = "indigo",linewidth=3)

plt.legend()     


     
#plt.savefig(savepath + 'displacement_fcm5lc203_cc.png',dpi=300,bbox_inches='tight')

plt.show()  

# cross section plot
for n in range(0,len(r_init)):
    if ( (long_init[n] < -112) & (long_init[n] > -112.5) ):
        plt.arrow(lat_init[n],r_init[n],lat_diff[n],r_diff[n],width=0.04,head_length=500)
        
plt.ylabel('depth (m)')    
#plt.savefig(savepath + 'displacementslice.png',dpi=300,bbox_inches='tight')   

# add colormap




# all horizontal in a specific depth range

# for n in range(0,len(r_init)):
#     if ( (r_init[n] < rmax) & (r_init[n] > rmin) ):
#         #plt.arrow(lat_init[n],r_init[n],lat_diff[n],r_diff[n],width=0.04,head_length=500)
#         plt.scatter(lat_init[n],r_init[n])
        
# plt.ylabel('depth (m)') 
# plt.title("All particles starting between depth "+str(rmin/1000)+"-"+str(rmax/1000)+"km")   
  
        
   
        
        