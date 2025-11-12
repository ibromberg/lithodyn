#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 16:04:18 2025

@author: ibromberg
"""

import numpy as np
from scipy.linalg import eig
from scipy.interpolate import griddata
import matplotlib.pyplot as plt

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

    return sigma_1, sigma_2, dir_1, dir_2

# ---------------

rEarth = 6371e3
sc = 4

# Read data
dat = np.loadtxt('stress_17ma_UCLC.txt', skiprows=9)

x = dat[:,0]
y = dat[:,1]
z = dat[:,2]
tau_cart = dat[:,3:12]  # 9 stress components

# loop through all time steps

# convert to spherical
long, lat, r = cart2sph(x, y, z)
long = -long
elev_km = (r - rEarth)/1000

# grid
minLong, maxLong = -126, -101
minLat, maxLat = 26, 44
ddeg = 0.5
mlong, mlat = np.meshgrid(np.arange(minLong, maxLong+ddeg, ddeg),
                          np.arange(minLat, maxLat+ddeg, ddeg))

theta = long
phi = lat

# principal stresses
sigma_1, sigma_2, dir_1, dir_2 = tau_cart2sphr(tau_cart, r, theta, phi)
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

# --------- PLOTTING

plt.figure()
plt.quiver(comp_long, comp_lat, comp_east, comp_north, sc, color='k')
#plt.quiver(tens_long, tens_lat, tens_east, tens_north, sc, color='b')
#plt.legend(['Compressive principle stress','Extensive'])
plt.show()