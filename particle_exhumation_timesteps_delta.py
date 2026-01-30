#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 16:00:01 2025

@author: ibromberg

exhumation rather than particle displacement

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
from scipy.ndimage import gaussian_filter

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

def grid(lat, long, r, lati, longi, sig):
    Long, Lat = np.meshgrid(longi, lati)

    R = griddata( (long, lat), r, (Long, Lat), method='linear')

    R_near = griddata((long, lat), r, (Long, Lat), method='nearest')
    R = np.where(np.isnan(R), R_near, R)
    R = gaussian_filter(R, sigma=sig)
    
    return Lat, Long, R

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
cclong, cclat, ccdeg = core['long'], core['lat'], core['deg']

# ----------------

varyLC = True

if varyLC == True:
    path = '/Users/ibromberg/Documents/COMSOL 63/agu2025/modelout/particles.txt'
    topopath = '/Users/ibromberg/Documents/COMSOL 63/agu2025/modelout/topo.txt'
    savepath = '/Users/ibromberg/Documents/COMSOL 63/prelim/varyLC/'
else:
    path = '/Users/ibromberg/Documents/COMSOL 63/particle data/fcm5lc21_allparticles_3ma.txt'
    topopath = '/Users/ibromberg/Documents/COMSOL 63/topodata/fcm5lc21_topo_all_3ma.txt'
    savepath = '/Users/ibromberg/Documents/COMSOL 63/prelim/lc21/'


#savepath = '/Users/ibromberg/Documents/COMSOL 63/plots/particles_anim/exhumation/'
#savepath = '/Users/ibromberg/Documents/COMSOL 63/agu2025/plots/exhumation_2/'


zoomed = True # make 'True' if you want it zoomed in on the core complexes
logscale = True # make 'True' to put the color bar on a log scale
plotalpha = False # true to make two plots for debugging - plots the change between each frame
stationary = True # if generating specific frames for plots in papers/writeups

ma_init = 17.00
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

rmin, rmax = -20000, 0
longmin, longmax = -122, -106
latmin, latmax = 28, 42

ccsize, ccline = 25, 3
cccolor = 'gold'

sig = 1 #smoothing when zoomed out

if zoomed == True:
    # core complexs
    longmin, longmax = -118, -110
    latmin, latmax = 33, 40
    ccsize, ccline = 85, 3
    sig = 2
    
n = 100 # grid topo    initial and final values

#when zoomed in, sig = 2; corresponds to smoothing over ~2 grid cells to suppress particle-scale noise; 
# smoothing over ~30–45 km removes particle discreteness & preserves regional gradients and structure

longi = np.linspace(longmin, longmax, n)
lati  = np.linspace(latmin,  latmax,  n)

if stationary == True: 
    plt.rcParams.update({'font.size': 28})
    contourfont = 24
    ccsize = 200
else:    
    plt.rcParams.update({'font.size': 18})
    contourfont = 16


# calculate the initial difference

df_init = df_particle[0].copy() 
df_topo = df_topos[0].copy()

long_init, lat_init, r_init = process(df_init)
long_t, lat_t, r_t = process(df_topo)

if zoomed == True:
    long_t, lat_t, r_t = zip(*[
        (x, y, z) for x, y, z in zip(long_t, lat_t, r_t)
        if x <= longmax and x >= longmin and y <= latmax and y >= latmin            
        ])
    
df_init["lat"] = lat_init    
df_init["long"] = long_init
df_init["r"] = r_init

df_init = df_init[(df_init.lat < latmax) & (df_init.long < longmax) & (df_init.long > longmin) & (df_init.r > rmin) & (df_init.r < rmax)]

lat_init, long_init, r_init = (np.array( df_init["lat"].values.tolist() ),
                                np.array( df_init["long"].values.tolist() ),
                                np.array( df_init["r"].values.tolist() ))

topo_lat0, topo_long0, topo_r0 = grid(lat_t, long_t, r_t, lati, longi, sig)
part_lat0, part_long0, part_r0 = grid(lat_init, long_init, r_init, lati, longi, sig)
exhume_0 = topo_r0 - part_r0 # this is alpha_0


# LOOP

fin = len(df_particle)-1

depth = 0

if stationary == True: start, step = 99, 100
else: start, step = 9, 10

for i in [9,99,199,299]:#range(start,fin,step): # len(df_particle) #start fin step
     
    df_start = df_particle[i].copy()
    df_topo_start = df_topos[i].copy()
    
    df_fin = df_particle[i+1].copy() # next time step       
    df_topo_fin = df_topos[i+1].copy()
    
    long_init, lat_init, r_init = process(df_start)
    long_t, lat_t, r_t = process(df_topo_start)
    
    long_fin, lat_fin, r_fin = process(df_fin)
    long_t_f, lat_t_f, r_t_f = process(df_topo_fin)
      
    #clip data if zoomed in for topography
    if zoomed == True: 
        long_t, lat_t, r_t = zip(*[
            (x, y, z) for x, y, z in zip(long_t, lat_t, r_t)
            if x <= longmax and x >= longmin and y <= latmax and y >= latmin            
            ])
        
        long_t_f, lat_t_f, r_t_f = zip(*[
            (x, y, z) for x, y, z in zip(long_t_f, lat_t_f, r_t_f)
            if x <= longmax and x >= longmin and y <= latmax and y >= latmin            
            ])
    
    # filter particle data
    # put particle data all back into dataframe for filtering
    df_start["lat"], df_start["long"], df_start["r"] = lat_init, long_init, r_init 
    df_fin["lat"], df_fin["long"], df_fin["r"] = lat_fin, long_fin, r_fin    
        
    df_start  = df_start[(df_start.lat < latmax) & (df_start.lat > latmin) & (df_start.long < longmax) & (df_start.long > longmin) & (df_start.r > rmin) & (df_start.r < rmax)]
    df_fin  = df_fin[(df_fin.lat < latmax) & (df_fin.lat > latmin) & (df_fin.long < longmax) & (df_fin.long > longmin) & (df_fin.r > rmin) & (df_fin.r < rmax)]
        
    #re-extract those filtered lists
    lat_init, long_init, r_init = df_start["lat"].to_numpy(), df_start["long"].to_numpy(), df_start["r"].to_numpy()
    lat_fin, long_fin, r_fin = df_fin["lat"].to_numpy(), df_fin["long"].to_numpy(), df_fin["r"].to_numpy()
    
    medr = np.median(r_init)
    # print(medr)
      
    topo_lat_init, topo_long_init, topo_r_init = grid(lat_t, long_t, r_t, lati, longi, sig)
    part_lat_init, part_long_init, part_r_init = grid(lat_init, long_init, r_init, lati, longi, sig)
    
    topo_lat, topo_long, topo_r = grid(lat_t_f, long_t_f, r_t_f, lati, longi, sig)
    part_lat, part_long, part_r = grid(lat_fin, long_fin, r_fin, lati, longi, sig)
    
    exhume_i = topo_r_init - part_r_init
    exhume_f = topo_r - part_r
    
    alpha = exhume_f - exhume_i
    exhume = exhume_i - exhume_0
    
    # find min and max for color scale checking
    print(np.max(exhume))
    print(np.min(exhume))
    
    if depth==0: depth = abs(medr)/1000
    
    
     # -------------------------------------------------------
    
     # plotting
     
    if plotalpha == True:
        fig, ax = plt.subplots(1, 2, dpi=300)
        fig.set_size_inches(20, 7)
        
        # deformed state lines
        a = 0.3
        lines = 'dashed'
        
        ax[0].plot(arizona[0],arizona[1], color = "black", ls = lines, alpha = a)
        ax[0].plot(california[0],california[1], color = "black", ls = lines, alpha = a)
        ax[0].plot(colorado[0],colorado[1], color = "black", ls = lines, alpha = a)
        #ax[0].plot(mexico[0],mexico[1], color = "black")
        ax[0].plot(nevada[0],nevada[1], color = "black", ls = lines, alpha = a)
        ax[0].plot(new_mexico[0],new_mexico[1], color = "black", ls = lines, alpha = a)
        ax[0].plot(oregon[0],oregon[1], color = "black", ls = lines, alpha = a)
        ax[0].plot(texas[0],texas[1], color = "black", ls = lines, alpha = a)
        ax[0].plot(utah[0],utah[1], color = "black", ls = lines, alpha = a)
        ax[0].plot(wyoming[0],wyoming[1], color = "black", ls = lines, alpha = a)
        
        ax[1].plot(arizona[0],arizona[1], color = "black", ls = lines, alpha = a)
        ax[1].plot(california[0],california[1], color = "black", ls = lines, alpha = a)
        ax[1].plot(colorado[0],colorado[1], color = "black", ls = lines, alpha = a)
        #ax[1].plot(mexico[0],mexico[1], color = "black")
        ax[1].plot(nevada[0],nevada[1], color = "black", ls = lines, alpha = a)
        ax[1].plot(new_mexico[0],new_mexico[1], color = "black", ls = lines, alpha = a)
        ax[1].plot(oregon[0],oregon[1], color = "black", ls = lines, alpha = a)
        ax[1].plot(texas[0],texas[1], color = "black", ls = lines, alpha = a)
        ax[1].plot(utah[0],utah[1], color = "black", ls = lines, alpha = a)
        ax[1].plot(wyoming[0],wyoming[1], color = "black", ls = lines, alpha = a)
        
        ax[0].set_xlim(longmin,longmax)
        ax[0].set_ylim(latmin,latmax)
        ax[0].set_xlabel("Longitude")
        ax[0].set_ylabel("Latitude")
        
        ax[1].set_xlim(longmin,longmax)
        ax[1].set_ylim(latmin,latmax)
        ax[1].set_xlabel("Longitude")
        ax[1].set_ylabel("Latitude")
        
        # plot topo minus particle     TOTAL    
        normalcolor1 = colors.TwoSlopeNorm(vmin=-5000,vcenter=0,vmax=5000) # for linear scale
        #normalcolor = colors.SymLogNorm(linthresh=10, linscale=1, vmin=-5000, vmax=5000, base=10) # for log scale
            
        ax[0].contourf(part_long,part_lat,exhume,cmap="seismic_r",levels=30, norm=normalcolor1) #,norm=normalcolor
        
        sm = ScalarMappable(norm=normalcolor1, cmap='seismic_r')
        sm.set_array([])  # Dummy array to satisfy matplotlib
        cbar=plt.colorbar(sm,ax=ax[0],label="Δ(Elevation - Particle Depth) (m)")
        
        # plot topo minus particle  ALPHA       
        normalcolor2 = colors.TwoSlopeNorm(vmin=-500,vcenter=0,vmax=500) # for linear scale
        #normalcolor = colors.SymLogNorm(linthresh=10, linscale=1, vmin=-5000, vmax=5000, base=10) # for log scale
            
        #alpha_clipped = np.clip(alpha, -500, 500)
        
        cf2=ax[1].contourf(part_long,part_lat,alpha,cmap="seismic",levels=30,norm=normalcolor2) #,norm=normalcolor
        #cbar=plt.colorbar(cf2,ax=ax[1],label="m")
        sm = ScalarMappable(norm=normalcolor2, cmap='seismic')
        sm.set_array([])  # Dummy array to satisfy matplotlib
        cbar=plt.colorbar(sm,ax=ax[1],label="Δ(Elevation - Particle Depth) (m)")
        
              
        # contour lines
        lev = [500, 1000, 1500, 2000, 2500, 3000, 3500]
        CS = ax[0].contour(topo_long, topo_lat, topo_r, levels=lev, colors='k',alpha=0.7)
        ax[0].clabel(CS, inline=True,fontsize=16,fmt='%d m')
        CS = ax[1].contour(topo_long, topo_lat, topo_r, levels=lev, colors='k',alpha=0.7)
        ax[1].clabel(CS, inline=True,fontsize=16,fmt='%d m')
        
        alphamax = alpha.max()
        alphamin = alpha.min()
        
        ma_diff = ma_init - (i+1)*timestep
        ma_diff_2 = ma_diff - timestep
        ax[0].set_title("Change in Elevation Minus Particle Depth from 17Ma to " + f"{ma_diff:.2f}" + "Ma, \n Average Particle Depth at 17Ma -" + f"{depth:.2f}" + "km ")
        ax[1].set_title("Change in Elevation Minus Particle Depth from " + f"{ma_diff:.2f}" + "Ma to " + f"{ma_diff_2:.2f}" + f"\n max: {alphamax:.2f}" + f", min: {alphamin:.2f}")
        
        #plt.savefig(savepath + str(i) + '.png',dpi=300) #,bbox_inches='tight'
        plt.show()
        
        #plt.scatter(part_lat,part_r,c=part_long)
        #plt.colorbar()
        
    else:     
        fig, ax = plt.subplots()
        fig.set_size_inches(12, 11)
        
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
        
        
        # plot topo minus particle        
        normalcolor = colors.TwoSlopeNorm(vmin=-5000,vcenter=0,vmax=5000) # for linear scale
        if logscale == True: normalcolor = colors.SymLogNorm(linthresh=100, linscale=1, vmin=-5000, vmax=5000, base=10) # for log scale
            
        plt.contourf(part_long,part_lat,exhume,cmap="RdBu",levels=100, norm=normalcolor) #,norm=normalcolor
        
        sm = ScalarMappable(norm=normalcolor, cmap='RdBu')
        sm.set_array([])  # Dummy array to satisfy matplotlib
        cbar=plt.colorbar(sm,ax=ax,label="Δ(Elevation - Particle Depth) (m)")
        
        #cbar=plt.colorbar(label="Vertical Particle Displacement (m)")
        
        # contour lines
        lev = [500, 1000, 1500, 2000, 2500, 3000, 3500]
        CS = plt.contour(topo_long, topo_lat, topo_r, levels=lev, colors='k',alpha=0.7)
        plt.clabel(CS, inline=True,fontsize=contourfont,fmt='%d m',colors='dimgray')
        
        ma_diff = ma_init - (i+1)*timestep
        #plt.title("Change in Elevation Minus Particle Depth from 17Ma to " + f"{ma_diff:.2f}" + "Ma, \nAverage Particle Depth at 17Ma -" + f"{depth:.2f}" + "km ")
        if stationary==True: plt.title("17Ma to " + f"{ma_diff:.2f}" + "Ma")
        else:
            plt.title("Change in Elevation Minus Particle Depth from \n17Ma to " + f"{ma_diff:.2f}" + "Ma\n")

        #print(ma_diff, exhume.max(), exhume.min())
        
         # core complexes
         # zoomed in: s=75, linewidth=3
         # zoomed out: s = 25, linewidth = 5
         
        # for j in range(0,len(ccdeg)):
        #     #draw_line(cclong[i],cclat[i],ccdeg[i],10)
        #     long_end, lat_end = draw_line(cclat[j],cclong[j],ccdeg[j],1)
        #     plt.plot([cclong[j], long_end],[cclat[j],lat_end],color = 'k',linewidth=ccline+1.5)
        #     plt.plot([cclong[j], long_end],[cclat[j],lat_end],color = cccolor,linewidth=ccline)
        
        plotcc = plt.scatter(cclong,cclat,marker='s',s=ccsize,color=cccolor,label="Core Complexes",edgecolors='k',zorder=10)
            
        plt.legend(loc = 'lower left')     
        
        #plt.savefig(savepath + str(i) + '.png',dpi=300,bbox_inches='tight') #,bbox_inches='tight'
            
        plt.show()
        plt.close()
        
        
        