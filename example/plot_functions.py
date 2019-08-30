# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 09:59:18 2019

@author: Alex Sternagel
"""

## Contains plot functions

import matplotlib.pyplot as plt
import numpy as np

## plot function for soil moisture profile
def plot_theta(theta, theta_init, theta_event, dim, z, time):
    
    plt.figure()

    # soil mositure profile, separated in the proportions of event and pre-event water    
    plt.plot(theta[0:dim-1,0] + theta_event[0:dim-1,0], z[0:dim-1,], 'r-', linewidth=2)
    plt.plot(theta[0:dim-1,0], z[0:dim-1], 'y-', linewidth=2) # only well-mixed water
    plt.plot(theta_init[0:dim-1,0], z[0:dim-1], 'k+',linewidth=2, markersize=5) # initial soil moisture profile
    plt.xlabel('theta [-]',fontsize=14)
    plt.ylabel('z [m]', fontsize=14)
    plt.title('Soil Water Content'  , fontsize=14)
    plt.legend(['Mixed+Event Water','Well-Mixed Water','Initial Moisture'], loc=4)
    plt.ylim(-1,0)
 
    return

## plot function for concentration profile
def plot_Cw(concfin, Cw, Cw_init, dim, theta, theta_init, z_plot):
    
    plt.figure()

    # soil mositure profile, separated in the proportions of event and pre-event water    
    plt.plot(Cw_init[0:dim-1,0]*theta_init[0:dim-1,0],z_plot[0:dim-1],'k-',linewidth=2)
    plt.plot((Cw[0:dim-1,0]*theta[0:dim-1,0]),z_plot[0:dim-1],'r-',linewidth=2)
    plt.plot(np.asarray(concfin)[:,0],-np.asarray(concfin)[:,2],'b-',linewidth=2)
    plt.xlabel('Concentration in Soil Phase [kg/m³]',fontsize=14)
    plt.ylabel('z [m]', fontsize=14)
    plt.title('Solute Concentration Profile at End of Simulation' , fontsize=14)
    plt.legend(['Initial Cw','Final Cw','Final Concentration Profile (Observed)','Final C-retdeg only in water','Final C-retdeg in entire layer'], loc=4)
    plt.ylim(-1,0)
    plt.xlim(0,0.03)
    
    return

## plot function for age profile
def plot_age(dim, mtx_avg_age, z_plot):
    
    plt.figure()
    
    plt.plot(mtx_avg_age[0:dim-1,0],z_plot[0:dim-1],'k-',linewidth=2)   
    plt.xlabel('Age [s]',fontsize=14)
    plt.ylabel('z [m]', fontsize=14)
    plt.title('Age Profile at End of Simulation' , fontsize=14)
    plt.legend(['Age'], loc=4)
    plt.ylim(-1,0)
    
    return

## plot function for mass balance and mass profile
def plot_mass(prec_int, prec_conc, concfin, theta, dim, Cw, dz, theta_init, Cw_init, z, z_plot):
    
    # Calculates mass balance
    mass_input = np.sum(prec_int * prec_conc * 600) * 1000 # total input mass
    mass_obs_final = np.sum(np.asarray(concfin)[:,0] * (np.asarray(concfin)[:,2] - np.asarray(concfin)[:,1]) * 1.4 * 1.4) * 1000 # final observed mass
    mass_sim_final = np.sum(theta[0:dim-1,0] * Cw[0:dim-1,0] * dz) * 1000 # final simulated mass

    recovery_sim = (mass_sim_final / mass_input) * 100 # recovery rate of simulated masses
    recovery_obs= (mass_obs_final / mass_input) * 100 # recovery rate of observed masses
    
    # mass plot
    plt.figure()
    
    plt.plot((theta_init[0:dim-1,0] * Cw_init[0:dim-1,0] * dz) * 1000,z_plot[0:dim-1],'k-',linewidth=2)
    plt.plot((theta[0:dim-1,0] * Cw[0:dim-1,0] * dz) * 1000,z_plot[0:dim-1],'r-',linewidth=2)
    plt.plot((np.asarray(concfin)[:,0] * (np.asarray(concfin)[:,2] - np.asarray(concfin)[:,1])) * 1000,-np.asarray(concfin)[:,2],'b-',linewidth=2)
    plt.xlabel('Mass [g]',fontsize=14)
    plt.ylabel('z [m]', fontsize=14)
    plt.title('Solute Mass Profile at End of Simulation' , fontsize=14)
    plt.legend(['Initial Mass Profile','Final Mass Profile (Simulated,conservative)','Final Mass Profile (Observed)','Final Mass-retdeg only in water','Final Mass-retdeg in entire layer'], loc=4)
    plt.ylim(-1,0)
    plt.xlim(0,6)
    
    return 