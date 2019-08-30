# -*- coding: utf-8 -*-
"""
Created on Fri Apr 12 13:46:35 2019

@author: Alex Sternagel
"""

import pandas as pd
import numpy as np

## Import own functions
from read_functions import *
from init_functions import *

## import data
# precipitation time and intensity     
prec_int, prec_time = read_precip (filepath="source/boundarycon/prec_int_specht.csv")

# precipitation concentration 
prec_conc = read_precip_conc (filepath="source/boundarycon/prec_conc_specht.csv")

# soil parameters
ks, ths, thr, alph, n_vg, stor, l_vg = read_soilpara (filepath="source/soils/soilpara_specht.dat")

# initial soil moisture state
theta_init = read_init_theta(filepath="source/initial_states/moist_init_specht.csv", z=z, dim=dim)

# initial concentration profile in matrix
Cw_init = read_init_Cw(filepath="source/initial_states/conc_init.csv", z=z, dim=dim)

# final observed tracer concentration profile of matrix
concfin = pd.read_csv("source/initial_states/conc_final_specht.csv", sep=' ', names=['Cw_final','start_depth','end_depth'])

## Preallocates arrays regarding the particles of both soil matrix and pfd and others

position_z_event = np.array([]) # position of event particles entering the soil matrix
age_event = np.array([np.nan]) # age of event particles particles entering the soil matrix
Cw_event = np.array([np.nan]) # concentration of event particles entering the soil matrix
t_mix_event = np.array([np.nan])

pfd_position_z = np.array([np.nan]) # position of particles within the pfd
pfd_age = np.array([np.nan]) # age of particles within the pfd
pfd_Cw_event = np.array([np.nan]) # concentration of particles entering pfd

m_surface = 0 # water surface storage
m_slt_surface = 0 # solute surface storage
m_input = 0

# just for testing mass consistency
mass_totalwaterinput = 0 # counter for amount of incoming precipitation water
mass_totalsoluteinput = 0 # counter for amount of incoming solute mass
    
bla = 0 # counter for water mass entering the soil matrix
bla1 = 0 # counter for water entering the pfd
   
solute_test = 0 # counter for solute mass entering the soil matrix
solute_test2 = 0 # counter for solute mass entering the pfd
    
theta_event = 0

## calculates initial states of psi and c using initial theta values
psi_init, c_init = psi_theta(alph=alph, dim=dim, n_vg=n_vg, stor=stor, theta=theta_init, thr=thr, ths=ths)

## calculates initial hydraulic conductivities using psi
k_init = k_psi(alph=alph, dim=dim, ks=ks, l_vg=l_vg, n_vg=n_vg, psi=psi_init)

## Creates lookup table for D and k values
prob = 1.0
D_table = np.zeros((1,nclass)) # lookup table for D
K_table = np.zeros((1,nclass)) # lookup table for k
theta_table = np.arange(thr,ths,((ths-thr)/nclass)).transpose() # array with soil moisture bins

#calculates for every soil moisture bin the matrix potential(psi_h), water capacity (c_h), diffusion coefficient (D_table) and hydraulic conductivity (K_table)
for i in range(0,nclass):
    theta_actual =  theta_table[i]* np.ones((dim,1))
    
    psi_h, c_h = psi_theta(alph=alph, dim=dim, n_vg=n_vg, stor=stor, theta=theta_actual, thr=thr, ths=ths)
    
    k_help = k_psi(alph=alph, dim=dim, ks=ks, l_vg=l_vg, n_vg=n_vg, psi=psi_h)
    
    if c_h[0] > 0:
       D_table[0][i] = k_help[0] / c_h[0]
    else: 
       D_table[0][i] = 0
       
    K_table[0][i] = k_help[0]

## Initialisation of particle tracking

# sets the working parameters to their initial values
k = k_init
c = c_init
psi = psi_init
theta = theta_init
Cw = Cw_init
pfd_Cw = pfd_Cw_initial

# initialises particle mass (m) distribution (n), positions (position_z) and concentrations (c_particle)
n, m = init_particles(theta=theta_init, dz=dz, dim=dim, N=N)

# initialises particle positions and initial particle concentrations, particle ages, retarded and degraded particle solute concentration
position_z, c_particle = init_particles_pos(z=z, dim=dim, n=n, Cw_init=Cw_init)
position_z = np.round(position_z,3)

age_particle = np.zeros((position_z.size,1))

position_z = pd.concat([pd.DataFrame(position_z), pd.DataFrame(age_particle), pd.DataFrame(c_particle)], axis=1)
position_z = np.asarray(position_z)

# time settings
time = prec_time[0] #start time (usually 0)
t_end = t_end + time
t_boundary = prec_time # time stepping of boundary data (time series of the precip_file)
i_time = 1

# initial input conditions at soil surface
Cw_eventparticles = prec_conc[0] # initial concentration of the new incoming particles
qb_u = -0.5 * (prec_int[0] + prec_int[1]) # flux density of precipitation water input

## Plot settings
z_plot = np.zeros(((z.size-1),1)) # Changes the increments of the simulated soil grid to those of the real observed one (for a better plot fitting)

for i in range (0, z.size-1):
    z_plot[i] = (z[i]+z[i+1]) / 2

z_plot = np.append(z[0], z_plot)
