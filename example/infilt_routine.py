# -*- coding: utf-8 -*-
"""
Created on Tue Jul  2 14:05:46 2019

@author: Alex Sternagel
"""

import numpy as np

## Import own functions
from model_functions import *

## Infiltration routine defining different scenarios for distribution of the precipitation water to the soil matrix and the pfd

if np.abs(qb_u) > 0 or m_surface > 0:

    ## Initialisation of input masses and surface storages

   m_surface = m_input + (np.abs(qb_u) * rho * 1 * dtc) + m_surface # rainwater initially stored in a virtual water surface storage
   m_input = 0
   m_slt_surface = m_slt_surface + (np.abs(qb_u) * 1 * dtc * Cw_eventparticles) # solute mass initially stored in a virtual solute surface storage
   Cw_eventparticles = m_slt_surface / (m_surface / rho)
 
   m_input1 = ((k[0] + ks) / 2) * ((np.abs(psi_init[0] - psi_init[1]) / dz[0]) + 1) * 1 * rho * dtc # water mass input into soil matrix
   m_input2 = (pfd_k * (np.pi * pfd_r**2)) * rho * dtc * n_mak  # water mass input into pfd

   mass_totalwaterinput = mass_totalwaterinput + (np.abs(qb_u) * rho * 1 * dtc) # control parameter, sums up the total incoming water mass
   mass_totalsoluteinput = mass_totalsoluteinput + (np.abs(qb_u) * dtc * Cw_eventparticles) # control parameter, sums up the total incoming solute mass

   if Cw_eventparticles == 0: # concentration input after rainfall event dependent on masses in surface storages
       Cw_eventparticles = m_slt_surface / (m_surface / rho) 
     

   theta_count = 0 # scenario counter 
   
   ## Scenario 1: Initially water and solutes flow into the soil matrix and residual overflow infiltrates simultaneously the pfd
 
   if np.round(theta[0],4) < (ths * 0.8):
     
     if m_input1 > m_surface:
        m_input1 = m_surface
     
     
     if Cw_eventparticles * (m_input1 / rho) > m_slt_surface:
        Cw_eventparticles = m_slt_surface / (m_input1 / rho)
     
     
     # calls function describing infiltration into soil matrix
     age_event, bla, Cw_event, m_input, m_slt_surface, m_surface, position_z_event, solute_test, theta, theta_count = infilt_mtx(t_mix_mean2=t_mix_mean2, t_mix_SD2=t_mix_SD2, t_mix_mean1=t_mix_mean1, t_mix_SD1=t_mix_SD1, t_mix_ratio=t_mix_ratio, age_event=age_event, bla=bla, Cw_event=Cw_event, Cw_eventparticles=Cw_eventparticles, dz=dz, m=m, m_input1=m_input1, m_slt_surface=m_slt_surface, m_surface=m_surface, position_z_event=position_z_event, rho=rho, solute_test=solute_test, theta=theta, time=time, z=z)
 

   # if the incoming precipitation mass is higher than the infiltrating mass(m_input1) it comes to a water accumulation of this residual water at the surface.This accumulation then flows immediately into the pfd with m_input2.
   if m_surface > 0:
     
     if m_input2 > m_surface:
        m_input2 = m_surface
     

     # calls function describing infiltration into pfd
     bla1, m_input, m_slt_surface, m_surface, pfd_position_z, pfd_theta, solute_test2, theta_count = infilt_pfd(bla1=bla1, Cw_eventparticles=Cw_eventparticles, m_input=m_input, m_input2=m_input2, m_slt_surface=m_slt_surface, m_surface=m_surface, n_mak=n_mak, pfd_dz=pfd_dz, pfd_m=pfd_m, pfd_position_z=pfd_position_z, pfd_r=pfd_r, pfd_theta=pfd_theta, pfd_z=pfd_z, rho=rho, solute_test2=solute_test2, theta_count=theta_count, time=time)
    
   ## Scenario 2: If the top grid element of the soil matris is saturated the water and solutes infiltrate only the pfd

   if theta_count == 0 and np.round(theta[0],2) == ths and pfd_theta[0] < 1:
    
      if m_input2 > m_surface:
         m_input2 = m_surface
    
      # calls function describing infiltration into pfd
      bla1, m_input, m_slt_surface, m_surface, pfd_position_z, pfd_theta, solute_test2, theta_count = infilt_pfd(bla1=bla1, Cw_eventparticles=Cw_eventparticles, m_input=m_input, m_input2=m_input2, m_slt_surface=m_slt_surface, m_surface=m_surface, n_mak=n_mak, pfd_dz=pfd_dz, pfd_m=pfd_m, pfd_position_z=pfd_position_z, pfd_r=pfd_r, pfd_theta=pfd_theta, pfd_z=pfd_z, rho=rho, solute_test2=solute_test2, theta_count=theta_count, time=time)
    
        

 