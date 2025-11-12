#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 23 12:19:51 2023

@author: ibromberg

creating my own differential velocity magnitude file

theta is long, phi is lat

use datagrid.py to get the velocities along a regular grid

big change apr 2025: add color map for vertical velocity (average values at surf and 28bsl?)
"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.cm import ScalarMappable

def cart2sph(x,y,z): # function for converting cartesian coords to spherical
    theta = np.arctan2(y,x) * 180/np.pi # azimuth
    phi = np.arctan2(z,np.sqrt(x**2 + y**2)) * 180/np.pi # elevation
    r = np.sqrt(x**2 + y**2 + z**2)#-6371000
    return theta, phi, r


def cart2sphV(long,lat,vx,vy,vz): # convert velocities to spherical
    theta = long * np.pi/180
    phi = [90-(i * np.pi/180) for i in lat] 

    # vr = ( vx * np.sin(phi) * np.cos(theta) + vy * np.sin(phi) 
    #         * np.sin(theta) + vz * np.cos(phi) )
    vphi = ( vx * np.cos(phi) * np.cos(theta) + vy * np.cos(phi) 
            * np.sin(theta) - vz * np.sin(phi) )
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

# -------------------------------------------------------

# read in data from comsol: topography, velocity topo, and velocity 28km bsl
# paste in pathnames


#yr = "17"
fcm = "5"
lc = "21"

#topo_str = "/Users/ibromberg/Documents/COMSOL 61/extrusion/text_outputs/fcm"+fcm+"lc"+lc+"_topo_"+yr+"ma.txt"
#v_topo_str = "/Users/ibromberg/Documents/COMSOL 61/extrusion/text_outputs/fcm"+fcm+"lc"+lc+"_v_topo_"+yr+"ma.txt"
#v_28_str = "/Users/ibromberg/Documents/COMSOL 61/extrusion/text_outputs/fcm"+fcm+"lc"+lc+"_v_28bsl_"+yr+"ma.txt"
#v_normal = "/Users/ibromberg/Documents/COMSOL 61/extrusion/text_outputs/fcm"+fcm+"lc"+lc+"_v_normal_"+yr+"ma.txt"

# read in data from comsol: topography, velocity topo, and velocity 28kmbsl
timestep = "17_0" # file name
timestep_t = "17.00"
path = '/Users/ibromberg/Documents/COMSOL 63/velocitydata/fcm5lc21_'
savepath = '/Users/ibromberg/Documents/COMSOL 63/plots/'

#'fcm2_nobc_isostacytest_wallsnoslip_bottomconstraint/v28km'+timestep+'.txt'

vtopo = pd.read_csv(path+'vtopo'+timestep+'.txt',
                    sep=" ",comment='%',header=None,names=["x","y","z","vx","vy","vz"])
v28 = pd.read_csv(path+'v28km'+timestep+'.txt',
                  sep=" ",comment='%',header=None,names=["x","y","z","vx","vy","vz"])
vnorm = pd.read_csv(path+'vnorm'+timestep+'.txt',
                  sep=" ",comment='%',header=None,names=["x","y","z","vn"])

# paste pathname here
savepath = '/Users/ibromberg/Documents/COMSOL 63/plots/diffv/'

#topo = pd.read_csv(topo_str ,sep=" ",
#                   header=None,names=["x","y","z","elev"])
#vtopo = pd.read_csv(vtopo,
#                    sep=" ",comment='%',header=None,names=["x","y","z","vx","vy","vz"])
#v28 = pd.read_csv(v28,
#                  sep=" ",comment='%',header=None,names=["x","y","z","vx","vy","vz"])
#vnorm = pd.read_csv(v_normal,
#                  sep=" ",comment='%',header=None,names=["x","y","z","vn"])

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
#x, y, z = (np.array(topo["x"].values.tolist()), 
#           np.array(topo["y"].values.tolist()), 
#           np.array(topo["z"].values.tolist()))
#theta, phi, r = cart2sph(x,y,z)

xnorm, ynorm, znorm, v_norm = (np.array(vnorm["x"].values.tolist()), 
           np.array(vnorm["y"].values.tolist()), 
           np.array(vnorm["z"].values.tolist()),
           np.array(vnorm["vn"].values.tolist()) )


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

vLongdiff = np.array(diff(vLongsurf,vLong28))
vLatdiff = np.array(diff(vLatsurf,vLat28))

vLongdiff = np.array(diff(vLong28,vLongsurf))
vLatdiff = np.array(diff(vLat28,vLatsurf))

#vLongdiff = np.mean( np.array([ vLongsurf, vLong28]), axis=0 )
#vLatdiff = np.mean( np.array([ vLatsurf, vLat28]), axis=0 )

# core complexes
cclong = core['long']#.values.tolist()
cclat = core['lat']#.values.tolist()
ccdeg = core['deg']

# -------------------------------------------------------
# 1 mm/yr = 3.171e-11 m/s

f, ax = plt.subplots(1)
f.set_size_inches(12, 10)

scaling = 5e-9

#plt.tricontourf(thetasurf,phisurf,v_norm,levels=50,cmap="PRGn")
#plt.colorbar(label="mm/yr")

# surface and 28km bsl velocities
#plt.tricontourf(theta,phi,r,cmap="gist_earth") # topography
# add scale=1e-9 to get them both on the same plotting scale
Q2 = plt.quiver(theta28,phi28, vLong28, vLat28,color='royalblue',label="28 km BSL Velocity",scale=scaling) # surface arrows
Q1 = plt.quiver(thetasurf,phisurf, vLongsurf, vLatsurf,color='firebrick',label="Surface Velocity",scale=scaling) # surface arrows


plt.quiverkey(Q1,0.2,0.35,3.171e-10,label="10 mm/yr",coordinates='figure',color='black')
plt.quiverkey(Q1,0.2,0.3,3.171e-10/2,label="5 mm/yr",coordinates='figure',color='black')
plt.quiverkey(Q1,0.2,0.25,3.171e-10/10,label="1 mm/yr",coordinates='figure',color='black')

# Q1 = plt.quiver(thetasurf,phisurf, vLongsurf, vLatsurf,color='royalblue',label="Surface Velocity",scale=0.5e-8) # surface arrows
# Q2 = plt.quiver(theta28,phi28, vLong28, vLat28,color='firebrick',label="28 km BSL Velocity",scale=0.5e-8) # surfa,ce arrows

# plt.quiverkey(Q1,0.3, 0.3, 3.171e-10,label="10 mm/yr",coordinates='figure')
# plt.quiverkey(Q2,0.3, 0.25, 3.171e-10,label="10 mm/yr",coordinates='figure')

plt.xlim([-122.5,-105])
plt.ylim([28,42])
plt.title("Velocities: "+timestep_t+" Ma",fontsize=30)

statecolor='black'
a=0.5

# add the deformed state lines
plt.plot(arizona[0],arizona[1], color = statecolor, alpha=a)
plt.plot(california[0],california[1], color = statecolor, alpha=a)
plt.plot(colorado[0],colorado[1], color =  statecolor, alpha=a)
#plt.plot(mexico[0],mexico[1], color = "black")
plt.plot(nevada[0],nevada[1], color =  statecolor, alpha=a)
plt.plot(new_mexico[0],new_mexico[1], color =  statecolor, alpha=a)
plt.plot(oregon[0],oregon[1], color =  statecolor, alpha=a)
plt.plot(texas[0],texas[1], color =  statecolor, alpha=a)
plt.plot(utah[0],utah[1], color = statecolor, alpha=a)
plt.plot(wyoming[0],wyoming[1], color =  statecolor, alpha=a)

# plot core complexes

# plot core complexes
for i in range(0,len(ccdeg)):
    #draw_line(cclong[i],cclat[i],ccdeg[i],10)
    long_end, lat_end = draw_line(cclat[i],cclong[i],ccdeg[i],1)
    plt.plot([cclong[i], long_end],[cclat[i],lat_end],color = "k",linewidth=4.5, zorder=10)
    plt.plot([cclong[i], long_end],[cclat[i],lat_end],color = "magenta",linewidth=3, zorder=10)

# dots = plt.scatter(mcc[0],mcc[1],marker='o',s=150,color="grey",edgecolors="black",linewidth=1,zorder=100)
plotcc = plt.scatter(cclong,cclat,marker='s',s=100,color='magenta',edgecolors='k',label="Core Complexes",zorder=11)


plt.legend()
plt.savefig(savepath + "surfand28v" + timestep_t + ".png",dpi=300,bbox_inches='tight')
#plt.show()

#plt.xlim([-112,-106])
#plt.ylim([34,40])
#plt.savefig(savepath1,dpi=300,bbox_inches='tight')

# differential velocities ---------------------------------------------

f, ax = plt.subplots(1)
f.set_size_inches(12, 10)


# COMMENT BELOW OUT if you want different scaled color bars for each time step
normalcolor = colors.TwoSlopeNorm(vmin=-5,vcenter=0,vmax=5) 

plt.tricontourf(thetasurf,phisurf,v_norm,levels=50,cmap="RdBu_r",norm=normalcolor)
#plt.colorbar(label="Surface Vertical Velocity (mm/yr)")

sm = ScalarMappable(norm=normalcolor, cmap='RdBu_r')
sm.set_array([])  # Dummy array to satisfy matplotlib

cbar=plt.colorbar(sm,ax=ax,label="Surface Vertical Velocity (mm/yr)")

# 8e-9 works for scaling too
scaling = 2e-9

Q3 = plt.quiver(thetasurf,phisurf, vLongdiff, vLatdiff,color='black',scale=scaling,label="Differential Velocity",zorder=12) # surface arrows
plt.quiverkey(Q3,0.2,0.3,3.171e-10,label="10 mm/yr",coordinates='figure')
plt.quiverkey(Q3,0.2,0.25,3.171e-10/2,label="5 mm/yr",coordinates='figure')
plt.quiverkey(Q3,0.2,0.2,3.171e-10/10,label="1 mm/yr",coordinates='figure')

plt.xlim([-122.5,-105])
plt.ylim([28,42])
plt.title("Differential Velocities: "+timestep_t+" Ma",fontsize=30)

# add the deformed state lines

plt.plot(arizona[0],arizona[1], color = statecolor,alpha=a)
plt.plot(california[0],california[1], color = statecolor,alpha=a)
plt.plot(colorado[0],colorado[1], color = statecolor,alpha=a)
#plt.plot(mexico[0],mexico[1], color = "black")
plt.plot(nevada[0],nevada[1], color = statecolor,alpha=a)
plt.plot(new_mexico[0],new_mexico[1], color = statecolor,alpha=a)
plt.plot(oregon[0],oregon[1], color = statecolor,alpha=a)
plt.plot(texas[0],texas[1], color = statecolor,alpha=a)
plt.plot(utah[0],utah[1], color = statecolor,alpha=a)
plt.plot(wyoming[0],wyoming[1], color = statecolor,alpha=a)

# plot core complexes
for i in range(0,len(ccdeg)):
    #draw_line(cclong[i],cclat[i],ccdeg[i],10)
    long_end, lat_end = draw_line(cclat[i],cclong[i],ccdeg[i],1)
    plt.plot([cclong[i], long_end],[cclat[i],lat_end],color = "k",linewidth=4.5, zorder=10)
    plt.plot([cclong[i], long_end],[cclat[i],lat_end],color = "magenta",linewidth=3, zorder=10)

# dots = plt.scatter(mcc[0],mcc[1],marker='o',s=150,color="grey",edgecolors="black",linewidth=1,zorder=100)
plotcc = plt.scatter(cclong,cclat,marker='s',s=100,color='magenta',edgecolors='k',label="Core Complexes",zorder=11)

plt.legend()

plt.savefig(savepath + "diffv" + timestep_t + ".png",dpi=300,bbox_inches='tight')

#zoom in on CC:
#plt.xlim([-116,-112])
#plt.ylim([33,36])
