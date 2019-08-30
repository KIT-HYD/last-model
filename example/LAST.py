# -*- coding: utf-8 -*-
"""
Created on Mon Apr 15 13:12:13 2019

@author: Alex Sternagel
"""

import os
import numpy as np

## set path of working directory
dirpath= os.path.dirname(__file__)
os.chdir(dirpath)

## Import own functions
from model_functions import *
from plot_functions import *

## Initialisation: Predefines global, soil matrix and pfd parameters and loads the initialisation script

# global parameters
rho = 1000 # water density
t_end = 86400 * 1 # end time of simulation (default: 1 day = 86400 s)
dtc = 120 # simulation time step

# soil matrix parameters
z = np.arange(0,-1.6,-0.1) # size of total grid
dz = np.abs(np.diff(z)) # size of one grid element
dim = z.size
mob_fak = 0.1 # mobile particle fraction within soil matrix 
N = 1000000 # total number of particles within soil matrix
nclass = 800 # number of bins to subdivide the diffusivity among the different particles

# tmix distribution
t_mix_ratio = 1.0 # to switch off tmix distribution, set t_mix_ratio to 1.0. Then, t_mix_mean1 is the only mixing time
t_mix_mean1 = 120
t_mix_mean2 = 12000
t_mix_SD1 = 60
t_mix_SD2 = 1200

# loading script defining geometry and parameters of macropore domain
exec(open('prefflow_domain.py').read())

# loading initialisation script with further parameters and input data
exec(open('init_LAST.py').read())

## Start of main model process with the routines for infiltration, displacement and mixing

while time < t_end:
    time = time + dtc

    # Infiltration routine
    exec(open('infilt_routine.py').read())
    
    # Displacement routine
    for pc in range(0,2):
    
        # initialisation of parameters v and D
        v = 1 * k # advective velocity in each grid element of the soil matrix
        D = 1 * k / c # diffusivity in each grid element of the soil matrix

        ipres = theta < 1.1 * thr # finds all grid elements with a soil moisture near to thr (almost dry soil)
        if ipres.any() == True: # sets the velocity and diffusivity in this grid elements to 0-> no flux!
            v[ipres == True] = 0
            D[ipres == True] = 0
   
        # displacement of pre-event particles in soil marix
        Cw, mtx_avg_age, position_z, theta = displ_mtx_pre(D=D, dim=dim, dtc=(0.5*dtc), D_table=D_table, dz=dz, K_table=K_table, m=m, mob_fak=mob_fak, position_z=position_z, ths=ths, v=v, z=z)
        
        # displacement of event particles in soil matrix
        if position_z_event.size > 0:
            position_z_event, theta_event = displ_mtx_event(dim=dim, dtc=(0.5*dtc), D_table=D_table, K_table=K_table, m=m, dz=dz, position_z_event=position_z_event, prob=prob, ths=ths, z=z)
            
        # displacement of particles in pfd and drainage; only one times in the predictor step because it is just dependent on constant advective veloscity
        if pc == 0 and np.isnan(pfd_position_z).all() == False:
            pfd_particles, pfd_Cw, pfd_theta, pfd_position_z = displ_pfd(n_mak=n_mak, pfd_dim=pfd_dim, pfd_dz=pfd_dz, pfd_m=pfd_m, pfd_n=pfd_n, pfd_position_z=pfd_position_z, pfd_r=pfd_r, pfd_z=pfd_z)

        # update of parameters psi, c and k for the corrector step
        psi, c = psi_theta(alph=alph, dim=dim, n_vg=n_vg, stor=stor, theta=theta, thr=thr, ths=ths)
        k = k_psi(alph=alph, dim=dim, ks=ks, l_vg=l_vg, n_vg=n_vg, psi=psi)
        
        
    # Mixing of event particle with pre-event particles within the first grid elements of soil matrix
    if position_z_event.size > 0:
        age_event, Cw_event, position_z, position_z_event = mixing_mtx(age_event=age_event, Cw_event=Cw_event, position_z=position_z, position_z_event=position_z_event, time=time)

    val_pos = -np.sort(-position_z[:,0])
    idx = np.argsort(-position_z[:,0])    
    position_z = np.concatenate((np.reshape(val_pos, (val_pos[:,].size,1)),np.reshape(position_z[idx,1],(position_z[idx,1].size,1)),np.reshape(position_z[idx,2],(position_z[idx,2].size,1))), axis=1)               
    
    # Mixing between pfd and soil matrix
    if np.isnan(pfd_position_z).all() == False:
        pfd_position_z, position_z = mixing_pfd_mtx(dtc=dtc, k=k, ks=ks, m=m, mak_sml=mak_sml, mak_mid=mak_mid, n_mak=n_mak, particle_contact_grid=particle_contact_grid, pfd_dim=pfd_dim, pfd_dz=pfd_dz, pfd_m=pfd_m, pfd_n=pfd_n, pfd_particles=pfd_particles, pfd_position_z=pfd_position_z, pfd_r=pfd_r, pfd_theta=pfd_theta, pfd_z=pfd_z, position_z=position_z, psi=psi, rate_big=rate_big, rate_mid=rate_mid, rate_sml=rate_sml, rho=rho, theta=theta, ths=ths, z=z)
   
    # Updates input conditions at top of soil
    ipos = np.where(t_boundary <= time)[0][-1] # finds actual position in boundary conditions time series
    qb_u = -0.5 * (prec_int[ipos,] + prec_int[ipos+1,]) # updates flux density of precipitation water input
    Cw_eventparticles = prec_conc[ipos,] # updates the concentration of precipitation water input
    
## Plots the final states of soil moisture, concentration and mass profile in soil matrix

# soil moisture profile
plot_theta(theta=theta, theta_init=theta_init, theta_event=theta_event, dim=dim, z=z, time=time)  

# concentration profile
plot_Cw(concfin=concfin, Cw=Cw, Cw_init=Cw_init, dim=dim, theta=theta, theta_init=theta_init, z_plot=z_plot)

# age profile
plot_age(dim=dim, mtx_avg_age=mtx_avg_age, z_plot=z_plot)

# mass balance and mass profile
plot_mass(prec_int=prec_int, prec_conc=prec_conc, concfin=concfin, theta=theta, dim=dim, Cw=Cw, dz=dz, theta_init=theta_init, Cw_init=Cw_init, z=z, z_plot=z_plot)





 
        
        
        
  

    