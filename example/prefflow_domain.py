# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 11:43:24 2019

@author: Alex Sternagel
"""

import pandas as pd
import numpy as np

#assumption: preferential flow domain with certain number of macropores shaped like a straight circular cylinder. 
#Water particles are spherical shaped in a cubic storage

## Initialisation of pfd and macropore geometry

pfd_z = np.arange(0,-1.05,-0.05) # maximum length of a macropore
pfd_D = (5/(10*100)) # diameter of macropore (mm)
n_mak = 16; # amount of macropores within the domain

pfd_r = pfd_D/2 # radius of a macropore
pfd_dz = np.abs(np.diff(pfd_z)) # length of a macropore grid element
pfd_dim = pfd_z.size
pfd_N = 10000 # amount of particles within pfd
rho = 1000 # density water 

pfd_theta = np.zeros((pfd_dim,1)) # initial soil moisture in pfd
pfd_Cw_initial = np.zeros((pfd_dim,1)) # initial solute concentration profile of pfd

## Initialisation of macropore depth distribution

mak_big = pfd_z[pfd_dim-1] # depth distribution of macropores within pfd
mak_mid = -0.8
mak_sml = -0.5

rate_big = 0.13 # distribution factors for the proportion of macropores and diffusive mixing masses
rate_mid = 0.19 # when adjusting the depth of the biggest macropores, please do it by adjusting "pfd_z" above
rate_sml = 0.68 

## Calculation of further pfd/macropore parameters

pfd_maxV_grid = np.pi * (pfd_r**2) * pfd_dz # volume of grid element
V_gesamt = pfd_maxV_grid.sum() # total volume of a macropore

M_eventgrid = rho * pfd_maxV_grid # water mass fitting in a grid element
m_gesamt = M_eventgrid.sum() # total water mass fitting in a macropore 
pfd_m = m_gesamt / pfd_N # mass of one particle

V_particle = pfd_m / rho # volume of a particle 
pfd_D_particle = (V_particle / (np.pi/6)) ** (1/3) # diameter of a particle
r_particle = pfd_D_particle / 2 # radius of a macropore

pfd_n = M_eventgrid / pfd_m # amount of particles fitting into a grid element
n_gesamt = pfd_n.sum() # amount of particles fitting into a macropore

U_eventgrid = 2 * np.pi * pfd_r # circumference of a grid element

particle_contact = np.floor(U_eventgrid/pfd_D_particle) # amount of particles fitting with their diameter next to each other in a row 
anzahl_reihen = np.floor(pfd_dz[0] / pfd_D_particle) # amount of rows fitting over each other 
particle_contact_grid = particle_contact * anzahl_reihen # amount of particles having contact to lateral surface of grid element
particle_contact_total = particle_contact_grid * pfd_dz.size # total amount of particles having contact to lateral surface of a macropore

if (particle_contact_total > pfd_N): # ensures that the amount of contact particles is not higher than the total amount of particles within the pfd
   particle_contact_total = pfd_N
   particle_contact_grid = particle_contact_total / pfd_dz.size

pfd_qmak = 2884.2 *(pfd_r**2) # flux density of a macropore
pfd_v = pfd_qmak # advective velocity of a particle 
pfd_k = pfd_v # hydraulic conductivity 
pfd_Diff = 0 # no diffusion within macropore
