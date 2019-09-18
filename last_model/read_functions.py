# -*- coding: utf-8 -*-
"""
Created on Wed Aug 28 10:02:42 2019

@author: Alex Sternagel
"""

## Contains read in-functions for input data

import pandas as pd
import numpy as np

## input of precipitation time series and intensity
def read_precip (prec):
    #prec = pd.read_csv(filepath, sep=' ', names=['prec_time','prec_int'])
    #prec = np.asarray(prec)

    prec_int = prec[:,1]
    prec_time = prec[:,0]

    conv_prec_int = (1/1000)/3600
    prec_int = prec_int * conv_prec_int

    conv_prec_time = 60
    prec_time = prec_time * conv_prec_time
    
    return prec_int, prec_time

## input of precipitation solute concentration
def read_precip_conc (prec_conc):
    #prec_conc = pd.read_csv(filepath, sep=' ', names=['prec_time','prec_conc'])
    #prec_conc = np.asarray(prec_conc)

    prec_conc = prec_conc[:,1]

    conv_conc = 1
    prec_conc = prec_conc * conv_conc
    
    return prec_conc

## input of soil parameters (only valid for one soil type)
def read_soilpara (filepath):
    raise DeprecationWarning('Please define the soil parameters in the config file')

    soil = pd.read_csv(filepath, sep=' ', names=['ks','ths','thr','alph','n_vg','stor', 'l_vg', 'depth1','depth2'])
    soil = np.asarray(soil)
    
    ks = soil[0][0]
    ths = soil[0][1]
    thr = soil[0][2]
    alph = soil[0][3]
    n_vg = soil[0][4]
    stor = soil[0][5]
    l_vg = soil[0][6]
    
    return ks, ths, thr, alph, n_vg, stor, l_vg

## input of initial soil moisture profile
def read_init_theta (theta, z, dim):
    #theta = pd.read_csv(filepath, sep=' ', names=['theta','start_depth','end_depth'])
    #theta = np.asarray(theta)
    m_i = theta.shape[0] # number or rows = number of observation values
    
    start_depth = np.zeros((m_i,1))
    end_depth = np.zeros((m_i,1))
    theta_init = theta[m_i-1][0] * np.ones((dim,1))
    
    start_depth[0] = theta[0][1]
    end_depth[0] = theta[0][2]
    
    
    for i in range(1,m_i):
        start_depth[i] = theta[i][1]
        end_depth[i] = theta[i][2]
        
        pos = np.where((z <= -start_depth[i-1]) & (z > -end_depth[i-1]))        
        pos = np.asarray(pos).transpose()
        
        for j in range (0,pos.shape[0]):
            theta_init[pos[j][0]] = theta[i-1][0] + (theta[i][0]-theta[i-1][0]) * ((z[pos[j][0]] * -1) - start_depth[i-1]) / (end_depth[i-1] - start_depth[i-1])
            
        
    return  theta_init

## input of initial concentration profile of matrix
def read_init_Cw (Cw, z, dim):
    #Cw = pd.read_csv(filepath, sep=' ', names=['Cw','start_depth','end_depth'])
    #Cw = np.asarray(Cw)
    
    m_i = Cw.shape[0] # number or rows = number of observation values
    
    start_depth = np.zeros((m_i,1))
    end_depth = np.zeros((m_i,1))
    Cw_init = Cw[m_i-1][0] * np.ones((dim,1))
    
    start_depth[0] = Cw[0][1]
    end_depth[0] = Cw[0][2]
    
    for i in range(1,m_i):
        start_depth[i] = Cw[i][1]
        end_depth[i] = Cw[i][2]
        
        pos = np.where((z <= -start_depth[i-1]) & (z > -end_depth[i-1]))        
        pos = np.asarray(pos).transpose()
        
        for j in range (0,pos.shape[0]):
            Cw_init[pos[j][0]] = Cw[i-1][0] + (Cw[i][0]-Cw[i-1][0]) * ((z[pos[j][0]] * -1) - start_depth[i-1]) / (end_depth[i-1] - start_depth[i-1])
        
    return  Cw_init