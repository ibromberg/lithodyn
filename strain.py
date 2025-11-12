#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 10 16:11:08 2025

@author: ibromberg

"""

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from scipy.interpolate import griddata
from scipy.ndimage import gaussian_filter

def cart2sph(x,y,z): # function for converting cartesian coords to spherical
    theta = np.arctan2(y,x) #* 180/np.pi# azimuth
    phi = np.arctan2(z,np.sqrt(x**2 + y**2)) #* 180/np.pi # elevation
    r = np.sqrt(x**2 + y**2 + z**2)#-6371000
    return r, theta, phi

def rot_sph(theta, phi):
    # rotation matrix A from Cartesian to spherical basis
    
    A = np.array([ [np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), np.cos(theta) ],
                   [np.cos(theta)*np.cos(phi), np.cos(theta)*np.sin(phi), -np.sin(theta)],
                   [-np.sin(phi),               np.cos(phi),               0            ]  ])
    return A

def strain_tensor_cart(e_xx, e_xy, e_xz, e_yx, e_yy, e_yz, e_zx, e_zy, e_zz):
   # make 3x3 cartesian strain rate tensor
    tensor =  np.array([ [e_xx, e_xy, e_xz],
                         [e_yx, e_yy, e_yz],
                         [e_zx, e_zy, e_zz]  ])
    
    return tensor

def rotate_strain_to_sph(E_cart, A):
    # rotate strain tensor from cart to sph: E_sph = A E A^T
    spherical = A @ E_cart @ A.T
    return spherical

R_earth = 6371000.0  # meters
# -----------

# columns: x, y, z, e_xx, e_xy, e_xz, e_yx, e_yy, e_yz, e_zx, e_zy, e_zz
# path = data file from comsol, title = for saving file, titlestring = for plot
path = '/Users/ibromberg/Documents/COMSOL 63/strain/strain_17_0.txt'
title = 'strainUCLC_17.png' 
titlestring = 'UCLC Boundary 17Ma' 

df = pd.read_csv(path, sep=" ", comment="%", header=None,
                 names=['x'  , 'y'  , 'z',
                        'exx', 'exy', 'exz',
                        'eyx', 'eyy', 'eyz',
                        'ezx', 'ezy', 'ezz']  )

# split all columns into arrays 
x, y, z = df[['x', 'y', 'z']].to_numpy().T
exx, exy, exz, eyx, eyy, eyz, ezx, ezy, ezz = df.iloc[:, 3:].to_numpy().T

# initialize spherical component lists

err, ert, erp = [], [], []
etr, ett, etp = [], [], []
epr, ept, epp = [], [], []

rp, thetap, phip = [], [], []

for i in range(0, len(df)):
    r, theta, phi = cart2sph(x[i],y[i],z[i])
    rp.append(r)
    thetap.append(theta)
    phip.append(phi)
    
    matrix_cart = strain_tensor_cart(exx[i], exy[i], exz[i],
                                     eyz[i], eyy[i], eyz[i],
                                     ezx[i], ezy[i], ezz[i]  )
    rot = rot_sph(theta, phi)
    spherical = rotate_strain_to_sph(matrix_cart, rot)
    
    err.append(spherical[0,0])
    ert.append(spherical[0,1])
    erp.append(spherical[0,2])
    
    etr.append(spherical[1,0])
    ett.append(spherical[1,1])
    etp.append(spherical[1,2])
    
    epr.append(spherical[2,0])
    ept.append(spherical[2,1])
    epp.append(spherical[2,2])
    
# err ert erp
# etr ett etp
# epr ept epp    looking at non-radial stress - ignore row 1 and column 1! 

# calculate  principal stresses: s1,2 = (epp+ett)/2 +- sqrt( (epp-ett)^2/2 + ept*etp )
#s1, s2 = [], []

#for i in range(0,len(ett)):
#    s1.append( (epp[i]+ett[i])/2 + np.sqrt( (epp[i]-ett[i])**2/2 + ept[i]*etp[i] ) )
#    s2.append( (epp[i]+ett[i])/2 - np.sqrt( (epp[i]-ett[i])**2/2 + ept[i]*etp[i] ) )

# doing this with x and y does fill the holes in but it also is meaningless data 
#for i in range(0,len(ett)):
#    s1.append( (exx[i]+eyy[i])/2 + np.sqrt( (exx[i]-eyy[i])**2/2 + exy[i]**2 ) )
#    s2.append( (exx[i]+eyy[i])/2 - np.sqrt( (exx[i]-eyy[i])**2/2 + exy[i]**2 ) )

if ept == etp:
    print('equal')
else:
    print('not equal')
    
r = np.array(rp)
ett = np.array(ett)
etp = np.array(etp)
ept = np.array(ept)
epp = np.array(epp)

#s1 = np.array(s1)
#s2 = np.array(s2)
theta = np.array(thetap)
phi = np.array(phip)

# make regular grid
n = 350 # the grid resolution
points = np.c_[theta, phi]

t_grid = np.linspace(min(theta), max(theta), n)
p_grid = np.linspace(min(phi),   max(phi),   n)

T, P = np.meshgrid(t_grid, p_grid)

# interpolate strain component
#Ftheta = griddata(points, s1, (T,P), method='linear')
#Fphi = griddata(points, s2, (T,P), method='linear')
Ett = griddata(points, ett, (T,P), method='linear')
Epp = griddata(points, epp, (T,P), method='linear')
Etp = griddata(points, etp, (T,P), method='linear')
Ept = griddata(points, ept, (T,P), method='linear')

#smooth for visualization
sigma = 2.0
Ett = gaussian_filter(Ett, sigma=sigma)
Epp = gaussian_filter(Epp, sigma=sigma)
Etp = gaussian_filter(Etp, sigma=sigma)
Ept = gaussian_filter(Ept, sigma=sigma)

# compute divergence using finite differences
dtheta = t_grid[1] - t_grid[0]
dphi = p_grid[1] - p_grid[0]

cosP = np.cos(P)
dEtt_dtheta = np.gradient(Ett *cosP, dtheta, axis=1)
dEpp_dphi   = np.gradient(Epp,       dphi  , axis=0)

div = (dEtt_dtheta + dEpp_dphi) / (R_earth * cosP)

#dFtheta_dt = np.gradient(Ftheta, dt, axis=1)
#dFphi_dp = np.gradient(Fphi, dp, axis=0)
#div = dFtheta_dt + dFphi_dp

# spherical ? looks worse, possibly calculated wrong
#div_sphere = (1 / np.sin(T)) * (dFtheta_dt * np.sin(T) + dFphi_dp)

# second invarient J2
J2 = np.sqrt(0.5 * ((Ett - Epp)**2 + (Etp)**2 + (Ept)**2))

#theta_f = T.ravel()
#phi_f = P.ravel()
#div_f = div.ravel()

#theta_f, phi_f, div_f = zip(*[
#(x, y, z) for x, y, z in zip(theta_f, phi_f, div_f) if not np.isnan(z)
#and x < -106 and x > -122# this prat removes bugged outliers
#])

# --------
# PLOTTING

# convert to radians for display
T_deg = np.degrees(T)
P_deg = np.degrees(P)

# adjust magnitude for s^-1 scale
vmin_div, vmax_div = -9e2, 9e2
vmin_j2, vmax_j2 = 0, 6e7


# divergence
plt.pcolormesh(T_deg,P_deg,div,cmap='seismic', vmin=vmin_div, vmax=vmax_div)
#plt.tricontourf(theta_f,phi_f,div_f,cmap='coolwarm')
plt.colorbar(label="Stress (N/m^2)")

plt.title(titlestring)
#plt.savefig('/Users/ibromberg/Documents/PLOTTING POSTPROCESSING/scrap/' + title,dpi=300,bbox_inches='tight')   




