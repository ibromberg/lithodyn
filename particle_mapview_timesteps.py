#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 16:09:20 2025

@author: ibromberg
"""
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib import colors
from matplotlib.cm import ScalarMappable

import matplotlib.tri as tri
from matplotlib.colors import SymLogNorm

from matplotlib import cm

from scipy.interpolate import griddata

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


#-----------

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

# ----------------

path = '/Users/ibromberg/Documents/COMSOL 63/particle data/fcm5lc21_allparticles.txt'
#path = '/Users/ibromberg/Documents/COMSOL 61/variable_LC/fcm3_lc30_particles.txt'
topopath = '/Users/ibromberg/Documents/COMSOL 63/topodata/fcm5lc21_topo_all.txt'

savepath = '/Users/ibromberg/Documents/COMSOL 63/plots/particles_anim/'
zoomed = True # make 'True' if you want it zoomed in on the core complexes
logscale = False # make 'True' to put the color bar on a log scale

ma_init = 17
timestep = 0.01
rounding = 2

# particle dataframe
df = pd.read_csv(path,sep=" ",comment='%',header=None)
df = df.drop(0, axis=1)
df_particle = [df.iloc[:, i:i+6] for i in range(0,df.shape[1], 6)]
# Now df_particle is a list of DataFrames:
# dfs[0] -> columns 0-2 
# dfs[1] -> columns 3-5
# dfs[2] -> columns 6-8
# dfs[3] -> columns 9-11 etc. 

for i in range(0,len(df_particle)):    
    df_particle[i].columns = ['x', 'y', 'z', 'vx', 'vy', 'vz']
    
    
# topography dataframe
df = pd.read_csv(topopath,sep=" ",comment="%",header=None)

# remove redundant xyz that happens when u extract from comsol
df = df.drop(0, axis=1)
df = df.drop(1, axis=1)
df = df.drop(2, axis=1)

df_topos = [df.iloc[:, i:i+3] for i in range(0,df.shape[1], 3)]   
for i in range(0,len(df_topos)):    
    df_topos[i].columns = ['x', 'y', 'z']
 
# ---------
# southwest: -110 to -118, below 30 to 37
# core complexes: -111,-116, 33 to 36

rmin = -20000
rmax = 0
longmin = -122
longmax = -106
latmin = 28
latmax = 42

ccsize = 25
ccline = 3
cccolor = 'magenta'


if zoomed == True:
    # core complexs
    longmin = -116
    longmax = -111
    latmin = 33
    latmax = 36
    ccsize = 75
    ccline = 3
    #cccolor = 'indigo'
    
for i in range(0,int((len(df_particle)-1))):
     
    df_init = df_particle[i].copy() # earlier time step
    df_fin = df_particle[i+1].copy() # next time step
    
    # topography at this time step
    df_topo = df_topos[i+1].copy()
    
    long_init, lat_init, r_init = process(df_init)
    long_fin, lat_fin, r_fin = process(df_fin)
    long_t, lat_t, r_t = process(df_topo)
    
    # get topo on grid for contour plot
    
    #clip data if zoomed in
    if zoomed == True:
        long_t, lat_t, r_t = zip(*[
            (x, y, z) for x, y, z in zip(long_t, lat_t, r_t)
            if x <= longmax and x >= longmin and y <= latmax and y >= latmin            
            ])
    
    
    n = 100 # grid
    longi = np.linspace(min(long_t), max(long_t), n)
    lati = np.linspace(min(lat_t), max(lat_t), n)
    Long, Lat = np.meshgrid(longi,lati)
    R = griddata((long_t,lat_t), r_t, (Long,Lat),method='cubic')
    
    
    # put particle data all back into dataframe for filtering
    df_init["lat"] = lat_init    
    df_init["long"] = long_init
    df_init["r"] = r_init
    df_fin["lat"] = lat_fin
    df_fin["long"] = long_fin
    df_fin["r"] = r_fin
    
    df_fin  = df_fin[(df_init.lat < latmax) & (df_init.long < longmax) & (df_init.long > longmin) & (df_init.r > rmin) & (df_init.r < rmax)]
    df_init = df_init[(df_init.lat < latmax) & (df_init.long < longmax) & (df_init.long > longmin) & (df_init.r > rmin) & (df_init.r < rmax)]
    
    #re-extract those filtered lists
    
    lat_init, long_init, r_init = (np.array( df_init["lat"].values.tolist() ),
                                    np.array( df_init["long"].values.tolist() ),
                                    np.array( df_init["r"].values.tolist() ))
    lat_fin, long_fin, r_fin = (np.array( df_fin["lat"].values.tolist() ),
                                np.array( df_fin["long"].values.tolist() ),
                                np.array( df_fin["r"].values.tolist() ))
    
    r_diff = (r_fin - r_init)
    lat_diff = (lat_fin-lat_init)
    long_diff = (long_fin-long_init)
    medr = np.median(r_init)
    # print(medr)
    
    # remove nans
    r_diff = r_diff.tolist()
    long_init = long_init.tolist()
    lat_init = lat_init.tolist()
    
    r_diff, lat_init, long_init, lat_diff, long_diff, lat_fin, long_fin = zip(*[
    (u, v, x, y, z, m, n) for u, v, x, y, z, m, n in zip(r_diff, lat_init, long_init, lat_diff, long_diff, lat_fin, long_fin) if not np.isnan(u)
    and u <= 5000 and u >= -5000# this prat removes bugged outliers
    ])
     
    print(i, max(r_diff))
    
    # clip data if zoomed in
    if zoomed == True:
        r_diff, lat_init, long_init, lat_diff, long_diff, lat_fin, long_fin = zip(*[
        (u, v, x, y, z, m, n) for u, v, x, y, z, m, n in zip(r_diff, lat_init, long_init, lat_diff, long_diff, lat_fin, long_fin)
        if v <= latmax and v >= latmin and x <= longmax and x >= longmin
        ])
        
     # -------------------------------------------------------
    
     # plotting
    fig, ax = plt.subplots()
    fig.set_size_inches(12, 10)
    
     # add the deformed state lines
    a = 0.3
    lines = 'dashed'
    
    plt.plot(arizona[0],arizona[1], color = "black", ls = lines, alpha = a)
    plt.plot(california[0],california[1], color = "black", ls = lines, alpha = a)
    plt.plot(colorado[0],colorado[1], color = "black", ls = lines, alpha = a)
    #plt.plot(mexico[0],mexico[1], color = "black")
    plt.plot(nevada[0],nevada[1], color = "black", ls = lines, alpha = a)
    plt.plot(new_mexico[0],new_mexico[1], color = "black", ls = lines, alpha = a)
    plt.plot(oregon[0],oregon[1], color = "black", ls = lines, alpha = a)
    plt.plot(texas[0],texas[1], color = "black", ls = lines, alpha = a)
    plt.plot(utah[0],utah[1], color = "black", ls = lines, alpha = a)
    plt.plot(wyoming[0],wyoming[1], color = "black", ls = lines, alpha = a)
    
    plt.xlim(longmin,longmax)
    plt.ylim(latmin,latmax)
    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    
    
     #plt.axvline(x = -112, color = 'black',linestyle='dashed')
     #plt.axvline(x = -112.5, color = 'black',linestyle='dashed')
     
     #print(min(r_diff))
     #print(max(r_diff))
     
     # set color bar 0 at 0
    if i==0:
        normalcolor = colors.TwoSlopeNorm(vmin=min(r_diff),vcenter=0,vmax=max(r_diff)) 
        
        scaling = 0.3
        
    else:
        #
        if logscale == True:
            normalcolor = colors.SymLogNorm(linthresh=1, linscale=1, vmin=-25, vmax=25, base=10) # for log scale
        else:
            normalcolor = colors.TwoSlopeNorm(vmin=-50,vcenter=0,vmax=50) # for linear scale
        scaling = 0.01 #5e-9 # change as needed; 0.01 
         
    # Triangulation
    #triang = tri.Triangulation(long_init, lat_init)

   # Create a symmetric log normalization
    #normalcolor = SymLogNorm(linthresh=1, linscale=1, vmin=r_diff.min(), vmax=r_diff.max(), base=10)
   
   
    #plt.scatter(long_init,lat_init,c=r_diff, ec='black', linewidth=0.5, cmap='bwr', s=200, norm=normalcolor)
    #print(np.isnan(np.min(r_diff)))
    
    
    plt.tricontourf(long_init,lat_init,r_diff,cmap="seismic",norm=normalcolor)
    
    # contour lines
    lev = [500, 1000, 1500, 2000, 2500, 3000]
    CS = plt.contour(Long, Lat, R, levels=lev, colors='k',alpha=0.7)
    plt.clabel(CS, inline=True,fontsize=16,fmt='%d m')
    
    # to-scale arrows
    #for n in range(0,len(r_diff)): 
    #    plt.arrow(long_init[n],lat_init[n],long_diff[n],lat_diff[n])
    
    # 1 mm/yr = 3.171e-11 m/s
    # 1 degree = 111,111 m; 1000 m = 111.111 deg
    
    plt.subplots_adjust(bottom=0.2)
    
    
    Q1 = plt.quiver(long_init,lat_init, long_diff, lat_diff,color='black',scale=scaling)
    
    #plt.quiverkey(Q1,0.2,0.23,3.171e-11,label="1 mm",coordinates='figure',color='black')
    if i==0:
        qk = ax.quiverkey(Q1, X=0.25, Y=-0.2, U=1/111111 * 1000, label='1000 m Horizontal Displacement',labelpos='E') # horizontal displacement key
    else:
        qk = ax.quiverkey(Q1, X=0.25, Y=-0.2, U=1/111111 * 50, label='50 m Horizontal Displacement',labelpos='E') # horizontal displacement key
    #qk.set_bbox(dict(facecolor='white', edgecolor='black', boxstyle='round,pad=0.3'))
    
    sm = ScalarMappable(norm=normalcolor, cmap='seismic')
    sm.set_array([])  # Dummy array to satisfy matplotlib
    
    cbar=plt.colorbar(sm,ax=ax,label="Vertical Particle Displacement (m)")
    #cbar=plt.colorbar(label="Vertical Particle Displacement (m)")
    
    if i==0: depth = abs(medr)/1000
     #print(depth)
     
     #plt.title("Particle Displacement from 16Ma - 17Ma, Median Particle Start Depth " + str(round(depth,1)) + "km")
    ma_diff = ma_init - timestep
    plt.title("Particle Displacement from " + f"{ma_init:.2f}" + " Ma to " + f"{ma_diff:.2f}" + "Ma, Average Particle Depth at 17Ma -" + f"{depth:.2f}" + "m ")
    
    
     # core complexes
     # zoomed in: s=75, linewidth=3
     # zoomed out: s = 25, linewidth = 5
     
    for j in range(0,len(ccdeg)):
        #draw_line(cclong[i],cclat[i],ccdeg[i],10)
        long_end, lat_end = draw_line(cclat[j],cclong[j],ccdeg[j],1)
        plt.plot([cclong[j], long_end],[cclat[j],lat_end],color = 'k',linewidth=ccline+1.5)
        plt.plot([cclong[j], long_end],[cclat[j],lat_end],color = cccolor,linewidth=ccline)
    
    plotcc = plt.scatter(cclong,cclat,marker='s',s=ccsize,color=cccolor,label="Core Complexes",edgecolors='k',zorder=10)
    
    
    plt.legend(loc = 'lower left')     
    
    plt.savefig(savepath + str(i) + " " + f"{ma_init:.2f}" + '.png',dpi=300) #,bbox_inches='tight'
    
    
    plt.show()
     
    ma_init = ma_diff
    
    
    
    
