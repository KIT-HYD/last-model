# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 10:06:36 2019

@author: Alex Sternagel
"""

## Contains initilisation functions 

import numpy as np

## calculation of initial psi and c values
def psi_theta (alph, dim, n_vg, stor, theta, thr, ths):
    psi = np.zeros((dim,1))
    c = np.zeros((dim,1))
    y = np.zeros((dim,1))
    
    for i in range(0,dim):
        m_vg = np.round(1 - (1 / n_vg),4)
        th_star = (theta[i,0] - thr) / (ths - thr)
        psi[i] = -1 * ((1 - (th_star**(1/m_vg))) / (th_star**(1/m_vg))) ** (1/n_vg) / alph
        
        if theta[i,0] <= 0.99 * ths:
            y[i] = (-m_vg * (1 / (1 + (psi[i] * alph * -1) ** n_vg)) ** (m_vg+1) * n_vg * (psi[i] * alph * -1) ** (n_vg-1) * alph)
            c[i] = (-(ths - thr) * y[i])
        elif theta[i,0] > 0.99 * ths:
            c[i] = stor
            
    return psi, c
            
## calculation of initial hydraulic conductivities
def k_psi (alph, dim, ks, l_vg, n_vg, psi):
    k = np.zeros((dim,1))
    
    for i in range(0,dim):
        m_vg = 1 - (1 / n_vg)
        v = 1 + (alph * (np.abs(psi[i]))) ** n_vg
        k[i] = (ks * (1 - (1 - (v**(-1))) ** m_vg) ** 2 / (v ** (m_vg * l_vg)))
        
    return k

## initialisation of matrix particle mass and distribution
def init_particles (theta, dz, dim, N):
    M = np.ones((dim-1,1))
    n = np.ones((dim-1,1))
    
    # Calculates particles mass based on soil water content in grid elements distributes the particles over the depth according to the soil moisture
    for i in range(0,(dim-1)):
        M[i] = theta[i] * dz[i] * 1 * 1000 # mass stored in each layer (kg)
        n[i] = np.round(theta[i] / np.sum(theta) * N) # particle number in each layer
    
    m = np.sum(M) / N # mass of one particle
    theta_p = n * m / (dz[i] * 1 * 1000)
    m = theta[0] / theta_p[0] * m
    
    return n, m
    
## initialisation of matrix particle positions and solute concentrations
def init_particles_pos(z, dim, n, Cw_init):
    position_z = np.ones((dim-1,1)).tolist()
    position_z_new = np.array([])
    
    c_particle = np.ones((dim-1,1)).tolist()
    c_particle_new = np.array([])
    
    # Calculates positions and concentrations of every particle
    for i in range(0,dim-1):
        increment = ((z[i] * -1) + z[i+1]) / n[i]
        position_z[i] = np.arange(z[i],(z[i+1]),increment) # normal distribution of particles in each grid element
        c_particle[i] = Cw_init[i] * np.ones((int(n[i]),1)) # every particle in the matrix gets the concentration of its grid element
    
    position_z = np.asarray(position_z)
    c_particle = np.asarray(c_particle)
    
    for i in range(0,dim-1):
        container = position_z[i]
        position_z_new = np.append(position_z_new, container)
        
        container = c_particle[i]
        c_particle_new = np.append(c_particle_new, container)
        
    
    return position_z_new, c_particle_new
