#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 16:04:18 2025

@author: ibromberg

converted from Sarah Bischoff/Chris Calvelage matlab code

takes in 9-component stress tensor and plots principle components
"""

import numpy as np
from scipy.linalg import eig
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter
import matplotlib.pyplot as plt
import pandas as pd

def cart2sph(x, y, z):
    # xyz to lat long r
    r = np.sqrt(x**2 + y**2 + z**2)
    phi = np.degrees(np.arccos(z / r))           # polar angle
    theta = np.degrees(np.arccos(x / np.sqrt(x**2 + y**2)))  # azimuth
    long = theta
    lat = 90 - phi
    return long, lat, r


def tau_cart2sphr(tau_cart, r, theta, phi):
    # stress from cartesian to spherical, find principle stresses & orientations

    n = len(r)
    sigma_1 = np.zeros(n)
    sigma_2 = np.zeros(n)
    dir_1 = np.zeros((n, 2))
    dir_2 = np.zeros((n, 2))

    tau_xx = tau_cart[:,0]
    tau_xy = tau_cart[:,1]
    tau_xz = tau_cart[:,2]
    tau_yx = tau_cart[:,3]
    tau_yy = tau_cart[:,4]
    tau_yz = tau_cart[:,5]
    tau_zx = tau_cart[:,6]
    tau_zy = tau_cart[:,7]
    tau_zz = tau_cart[:,8]

    for ii in range(n):
        # Transformation matrices
        sphi, cphi = np.sin(np.radians(phi[ii])), np.cos(np.radians(phi[ii]))
        stheta, ctheta = np.sin(np.radians(theta[ii])), np.cos(np.radians(theta[ii]))

        T1 = np.array([
            [sphi*ctheta, sphi*stheta, cphi],
            [cphi*ctheta, cphi*stheta, -sphi],
            [-stheta,     ctheta,       0]
        ])
        S_cart = np.array([
            [tau_xx[ii], tau_xy[ii], tau_xz[ii]],
            [tau_yx[ii], tau_yy[ii], tau_yz[ii]],
            [tau_zx[ii], tau_zy[ii], tau_zz[ii]]
        ])
        T2 = np.array([
            [sphi*ctheta, cphi*ctheta, -stheta],
            [sphi*stheta, cphi*stheta, ctheta],
            [cphi,       -sphi,       0]
        ])

        S_sph = T1 @ S_cart @ T2

        # Keep only theta-phi submatrix
        S = S_sph[1:,1:]

        # Principal stresses
        vals, vecs = eig(S)
        vals = np.real(vals)
        sigma_1[ii] = np.max(vals)
        sigma_2[ii] = np.min(vals)

        imax = np.argmax(vals)
        imin = np.argmin(vals)

        dir_1[ii,:] = [-vecs[0,imax], vecs[1,imax]]
        dir_2[ii,:] = [-vecs[0,imin], vecs[1,imin]]
        
        #if j == 1: print(S)
       # np.dot(eigvecs[:,0], eigvecs[:,1])  # should be ~0
       # np.dot(eigvecs[:,0], eigvecs[:,2])  # should be ~0
       # np.dot(eigvecs[:,1], eigvecs[:,2])  # should be ~0
        rotation = np.array(S_sph[0,1]) - np.array(S_sph[1,0])
        
        
    return sigma_1, sigma_2, dir_1, dir_2, rotation

def norm(comp_e, comp_n, tens_e, tens_n):
    mini = min(comp_e.min(), comp_n.min(), tens_e.min(), tens_n.min())
    maxi = max(comp_e.max(), comp_n.max(), tens_e.max(), tens_n.max())
    
    #mini = stacked.min()
    #maxi = stacked.max()
    
    comp_e_n = (comp_e - mini) / (maxi-mini)
    comp_n_n = (comp_n - mini) / (maxi-mini)
    tens_e_n = (tens_e - mini) / (maxi-mini)
    tens_n_n = (tens_n - mini) / (maxi-mini)
    
    return comp_e_n, comp_n_n, tens_e_n, tens_n_n
                        


def process(df):
    x = np.array(df["x"].values.tolist())
    y = np.array(df["y"].values.tolist())
    z = np.array(df["z"].values.tolist())
    
    theta, phi, r = cart2sph(x, y, z)
    
    
    return theta, phi, r

# ---------------

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

rEarth = 6371e3

# Read data

path = '/Users/ibromberg/Documents/COMSOL 63/strain/stress_all.txt'
path = '/Users/ibromberg/Documents/COMSOL 61/extrusion/text_outputs/fcm5lc21_10bsl_stress.txt'
savepath = '/Users/ibromberg/Documents/COMSOL 63/stressplots/'
#path = '/Users/ibromberg/Documents/COMSOL 63/strain/stress_17ma_UCLC.txt'
topopath = '/Users/ibromberg/Documents/COMSOL 63/topodata/fcm5lc21_topo_all.txt'

normalize = False
zoomed = False

# grid variables
minLong, maxLong = -126, -102 # -126, -101
minLat, maxLat = 28, 43 # 26, 44
ddeg = 0.5

if zoomed == True:
    # core complexs
    minLong, maxLong = -118, -110
    minLat, maxLat = 33, 40
    ddeg = 0.5

ma_init = 17.00

# read in dataframe
dfs = pd.read_csv(path,sep=" ",comment='%',header=None)

# remove redundant xyz that happens when u extract from comsol
dfs = dfs.drop(0, axis=1)
dfs = dfs.drop(1, axis=1)
dfs = dfs.drop(2, axis=1)

df_s = [dfs.iloc[:, i:i+12] for i in range(0,dfs.shape[1], 12)]

for i in range(0,len(df_s)):    
    df_s[i].columns = ['x'  , 'y'  , 'z',
                       'exx', 'exy', 'exz',
                       'eyx', 'eyy', 'eyz',
                       'ezx', 'ezy', 'ezz']
    
        
# topography dataframe
df = pd.read_csv(topopath,sep=" ",comment="%",header=None)

# remove redundant xyz that happens when u extract from comsol
df = df.drop(0, axis=1)
df = df.drop(1, axis=1)
df = df.drop(2, axis=1)

df_topos = [df.iloc[:, i:i+3] for i in range(0,df.shape[1], 3)]   
for i in range(0,len(df_topos)):    
    df_topos[i].columns = ['x', 'y', 'z']

fin = len(df_s)
for j in range(0,fin): #len(df_s)
    
    time = ma_init-(j*0.01)
    
    # topo
    df_topo = df_topos[j].copy()
    long_t, lat_t, r_t = process(df_topo)
    long_t = -long_t
    r_t = (r_t - rEarth)
    
    # get topo on grid for contour plot
    
    #clip data if zoomed in for topography
    if zoomed == True:
        long_t, lat_t, r_t = zip(*[(x, y, z) for x, y, z in zip(long_t, lat_t, r_t)
            if x <= maxLong and x >= minLong and y <= maxLat and y >= minLat            
            ])
    
    
    n = 100 # grid topo    
    sig = 2
    # topo grid
    longi = np.linspace(min(long_t), max(long_t), n)
    lati = np.linspace(min(lat_t), max(lat_t), n)
    Long, Lat = np.meshgrid(longi,lati)
    R = griddata((long_t,lat_t), r_t, (Long,Lat),method='nearest')
    R = gaussian_filter(R, sigma=sig)
    
    # stress
    # split all columns into arrays 
    df = df_s[j].copy()
    x, y, z = df[['x', 'y', 'z']].to_numpy().T
    exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz = df.iloc[:, 3:].to_numpy().T
    
    tau_cart_a = [exx, exy, exz,
                eyx, eyy, eyz,
                ezx, ezy, ezz]
    tau_cart = np.column_stack(tau_cart_a)
        
    # loop through all time steps
    
    # convert to spherical
    long, lat, r = cart2sph(x, y, z)
    long = -long
    elev_km = (r - rEarth)/1000
    
    #print(min(long),max(long))
    #print(min(lat),max(lat))
    
    # grid
    mlong, mlat = np.meshgrid(np.arange(minLong, maxLong+ddeg, ddeg),
                              np.arange(minLat, maxLat+ddeg, ddeg))
    
    theta = long
    phi = lat
    
    # principal stresses
    sigma_1, sigma_2, dir_1, dir_2, rot = tau_cart2sphr(tau_cart, r, theta, phi)
    sigma_1 /= 1e7
    sigma_2 /= 1e7
    
    # interpolation (scatteredInterpolant in matlab, griddata here)
    points = np.column_stack((long, lat))
    m_sigma_1 = griddata(points, sigma_1, (mlong, mlat), method='linear')
    m_dir_1north = griddata(points, dir_1[:,0], (mlong, mlat), method='linear')
    m_dir_1east = griddata(points, dir_1[:,1], (mlong, mlat), method='linear')
    m_sigma_2 = griddata(points, sigma_2, (mlong, mlat), method='linear')
    m_dir_2north = griddata(points, dir_2[:,0], (mlong, mlat), method='linear')
    m_dir_2east = griddata(points, dir_2[:,1], (mlong, mlat), method='linear')
    
    # separate compressive vs tensile stresses
    sigma = np.concatenate([m_sigma_1.ravel(), m_sigma_2.ravel()])
    dir_north = np.concatenate([m_dir_1north.ravel(), m_dir_2north.ravel()])
    dir_east = np.concatenate([m_dir_1east.ravel(), m_dir_2east.ravel()])
    long_all = np.concatenate([mlong.ravel(), mlong.ravel()])
    lat_all = np.concatenate([mlat.ravel(), mlat.ravel()])
    
    indC = sigma < 0
    indT = sigma > 0
    
    comp_east = sigma[indC] * dir_east[indC]
    comp_north = sigma[indC] * dir_north[indC]
    comp_long = long_all[indC]
    comp_lat = lat_all[indC]
    
    tens_east = sigma[indT] * dir_east[indT]
    tens_north = sigma[indT] * dir_north[indT]
    tens_long = long_all[indT]
    tens_lat = lat_all[indT]
    
    # double up stresses in opposite directions for plotting T shapes
    comp_east_opp = -1 * comp_east
    comp_north_opp = -1 * comp_north    
    tens_east_opp = -1 * tens_east
    tens_north_opp = -1 * tens_north
    
    comp_long_all = np.vstack((comp_long,comp_long))
    comp_lat_all  = np.vstack((comp_lat,comp_lat))
    comp_east_all = np.vstack((comp_east,comp_east_opp))
    comp_north_all = np.vstack((comp_north,comp_north_opp))
    
    tens_long_all = np.vstack((tens_long,tens_long))
    tens_lat_all  = np.vstack((tens_lat,tens_lat))
    tens_east_all = np.vstack((tens_east,tens_east_opp))
    tens_north_all = np.vstack((tens_north,tens_north_opp))
    
    
    # normalized
    
    comp_east_n, comp_north_n, tens_east_n, tens_north_n = norm(comp_east, comp_north, tens_east, tens_north)

    comp_east_opp_n = -1 * comp_east_n
    comp_north_opp_n = -1 * comp_north_n  
    tens_east_opp_n = -1 * tens_east_n
    tens_north_opp_n = -1 * tens_north_n
    
    comp_east_all_n = np.vstack((comp_east_n,comp_east_opp_n))
    comp_north_all_n = np.vstack((comp_north_n,comp_north_opp_n))
    
    tens_east_all_n = np.vstack((tens_east_n,tens_east_opp_n))
    tens_north_all_n = np.vstack((tens_north_n,tens_north_opp_n))
    
    # --------- PLOTTING
    
    fig, ax = plt.subplots()
    fig.set_size_inches(12, 10)
    
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

    plt.xlabel("Longitude")
    plt.ylabel("Latitude")
    plt.xlim(minLong,maxLong)
    plt.ylim(minLat,maxLat)
    
    
    plt.title("Stress at -10m " + f"{time:.2f}" + "Ma")
    
    
    
    
    
    # test arrows; just one in each direction
    #Q1 = plt.quiver(comp_long, comp_lat, comp_east, comp_north, color='r',scale=scaling)
    #Q2 = plt.quiver(tens_long, tens_lat, tens_east, tens_north, color='b',scale=scaling)
    
    #qk = ax.quiverkey(Q1, X=0.2, Y=0.1, U=5, label='5 Pa, Compressive principle stress',labelpos='N') # horizontal displacement key
    #qk = ax.quiverkey(Q2, X=0.2, Y=0.05, U=5, label='5 Pa, Extensive principle stress',labelpos='N') # horizontal displacement key
    
    #plt.legend(['Compressive principle stress','Extensive principle stress'])
    
    # contour lines
    lev = [1000, 2000, 3000, 4000]
    CS = plt.contour(Long, Lat, R, levels=lev, colors='k',alpha=0.7)
    plt.clabel(CS, inline=True,fontsize=16,fmt='%d m')
    
    
    if normalize == True:
        scaling = 40
        
        h = 5
        Q1 = plt.quiver(comp_long_all, comp_lat_all, comp_east_all_n, comp_north_all_n, color='r',scale=scaling,pivot='tip',headwidth=h)
        Q2 = plt.quiver(tens_long_all, tens_lat_all, tens_east_all_n, tens_north_all_n, color='b',scale=scaling, headwidth=h)
        
        #qk = ax.quiverkey(Q1, X=0.2, Y=0.1, U=5, label='5 Pa, Compressive principle stress',labelpos='N') # horizontal displacement key
        #qk = ax.quiverkey(Q2, X=0.2, Y=0.05, U=5, label='5 Pa, Extensive principle stress',labelpos='N') # horizontal displacement key
        
    else:
        
        scaling = 200
    
        Q1 = plt.quiver(comp_long_all, comp_lat_all, comp_east_all, comp_north_all, color='r',scale=scaling,pivot='tip')
        Q2 = plt.quiver(tens_long_all, tens_lat_all, tens_east_all, tens_north_all, color='b',scale=scaling)
        
        qk = ax.quiverkey(Q1, X=0.2, Y=0.1, U=5, label='5 Pa, Compressive principle stress',labelpos='N') # horizontal displacement key
        qk = ax.quiverkey(Q2, X=0.2, Y=0.05, U=5, label='5 Pa, Extensive principle stress',labelpos='N') # horizontal displacement key
       
        lev = [500, 1000, 1500, 2000, 2500, 3000, 3500]
        CS = plt.contour(Long, Lat, R, levels=lev, colors='k',alpha=1)
        plt.clabel(CS, inline=True,fontsize=16,fmt='%d m')
    
    plt.savefig(savepath + str(j) + " " + f"{time:.2f}" + '.png',dpi=300) #,bbox_inches='tight'
    
    plt.show()
    
    
    
    
    