# -*- coding: utf-8 -*-

## Contains crucial model functions (Infiltration, Displacement, Mixing)

import pandas as pd
import numpy as np
import copy

## matrix infiltration routine
def infilt_mtx(t_mix_mean2,t_mix_SD2,t_mix_mean1,t_mix_SD1,t_mix_ratio, age_event, bla, Cw_event, Cw_eventparticles, dz, m, m_input1, m_slt_surface, m_surface, position_z_event, rho, solute_test, theta, time, z):
    if m_input1 > m:
  
        precip_particles = np.floor(m_input1 / m).astype(int) # amount of event particles to be injected, converts the mass input to number of particles
        bla = bla + (precip_particles * m) # counter for amount of water totally infiltrating the soil matrix
        m_surface = m_surface - m_input1 # updates water surface storage
        m_slt_surface = m_slt_surface - (Cw_eventparticles * (m_input1 / rho)) # updates solute mass surface storage
        solute_test = solute_test + (Cw_eventparticles * (m_input1 / rho)) # counter for amount of solutes totally infiltrating the soil matrix within the first scenario
        m_input = m_input1 - (precip_particles * m) # calculates the remaining water particles, which where not considered by the "floor" function, and converts them back into a mass (in mm) to add it in the next time step. So no water input gets lost
        precip_position = z[0] * np.ones((precip_particles[0],1)) # in every time step the new infiltrating particles get the position z(1)
        theta_count = 1 # sets scenario counter to 1, thus the routine will not run through the subsequent scenarios
        theta[0] = np.round(theta[0] + (precip_particles * m) / (rho*dz[0]*1),4) # updates the soil moisture of the first grid element
        age_event = time * np.ones((precip_particles[0],1)) # matrix with the age of the infiltrating event particles, every new particle gets the value of its timestep when entering the domain
 

        t_mix_event = np.append(np.random.randn(np.floor(precip_particles[0] * t_mix_ratio).astype(int)) * t_mix_SD1 + t_mix_mean1, np.random.randn(np.ceil(precip_particles[0] * (1-t_mix_ratio)).astype(int)) * t_mix_SD2 + t_mix_mean2)
        t_mix_event = np.abs(t_mix_event[np.random.permutation(np.size(t_mix_event))].reshape((np.size(t_mix_event),1)))


        #adds the number of event particles entering the soil matrix with their concentration value to the Cw_event matrix
        if Cw_eventparticles == 0:
            Cw_event = np.zeros((precip_particles[0],1))     
        else: 
            Cw_event = Cw_eventparticles * np.ones((precip_particles[0],1))
     
     
        position_z_event_new = np.concatenate((precip_position, age_event, Cw_event, t_mix_event), axis=1)
        
        
        if position_z_event.size > 0:  
            position_z_event = np.append(position_z_event_new, position_z_event, axis=0)
        else:
            position_z_event = position_z_event_new

    return age_event, bla, Cw_event, m_input, m_slt_surface, m_surface, position_z_event, solute_test, theta, theta_count

## pfd infiltration routine
def infilt_pfd(bla1, Cw_eventparticles, m_input, m_input2, m_slt_surface, m_surface, n_mak, pfd_dz, pfd_m, pfd_position_z, pfd_r, pfd_theta, pfd_z, rho, solute_test2, theta_count, time):
    if (m_input2 > pfd_m) and (pfd_theta[0] < 1):
    
        pfd_theta_diff = 1 - pfd_theta[0] # detects the amount of water the first grid element of the pfd can gather in this time step
        precip_particles = np.floor(m_input2 / pfd_m).astype(int) # amount of event particles to be injected, converts the mass input to number of particles
        m_surface = m_surface - m_input2 # updates water surface storage
        m_input = m_input2 -(precip_particles * pfd_m) # calculates the remaining water particles, which where not considered by the "floor" function, and converts them back into a mass (in mm) to add it in the next time step. So no water input gets lost
        theta_count = 1 # sets scenario counter to 1, thus the routine will not run through the subsequent scenarios
        pfd_theta[0]=pfd_theta[0] + ((precip_particles * pfd_m) / ((rho * pfd_dz[0] * np.pi * (pfd_r**2)) * n_mak)) # updates the soil moisture of the first grid element

        # if the first grid element of the pfd is oversaturated due to a too high input mass (m_input2) then the surplus is considered byupdating the following parameters
        if pfd_theta[0] > 1:
           m_surface = m_surface + ((pfd_theta[0] - 1) * ((rho * pfd_dz[0] * np.pi * (pfd_r**2)) * n_mak))
           precip_particles = np.floor((pfd_theta_diff * ((rho * pfd_dz[0] * np.pi * (pfd_r**2)) * n_mak) / pfd_m)).astype(int)
           m_input2 = precip_particles * pfd_m
           pfd_theta[0] = 1
        
        # ensures that the infiltrating solute mass is not greater than the solute mass storage
        if Cw_eventparticles * (m_input2 / rho) > m_slt_surface:
           Cw_eventparticles = m_slt_surface / (m_input2 / rho) 
        
        m_slt_surface = m_slt_surface - (Cw_eventparticles * (m_input2 / rho)) # updates solute mass surface storage
        solute_test2 = solute_test2 + (Cw_eventparticles * (m_input2 / rho)) # counter for amount of solutes totally infiltrating the soil matrix within the first scenario

        bla1 = bla1 + (precip_particles * pfd_m) # counter for amount of water totally infiltrating the pfd
        pfd_precip_position = pfd_z[0] * np.ones((precip_particles[0],1)) # in every time step the new infiltrating particles get the position pfd_z(1)
        pfd_age = time * np.ones((precip_particles[0],1)) # matrix with the age of the infiltrating event particles, every new particle gets the value of its timestep when entering the domain

        #adds the number of event particles entering the pfd with their concentration value to the pfd_Cw_event matrix
        if Cw_eventparticles == 0:                    
            pfd_Cw_event = np.zeros((precip_particles[0],1))
        else:
            pfd_Cw_event = Cw_eventparticles * np.ones((precip_particles[0],1))
        
        pfd_position_z_new = np.concatenate((pfd_precip_position, pfd_age, pfd_Cw_event), axis=1)
        
        if np.isnan(pfd_position_z).any() == False:  
            pfd_position_z = np.append(pfd_position_z_new, pfd_position_z, axis=0)
        else:
            pfd_position_z = pfd_position_z_new
        
    return bla1, m_input, m_slt_surface, m_surface, pfd_position_z, pfd_theta, solute_test2, theta_count

## displacement routine for pre-event particles in matrix
def displ_mtx_pre(D, dim, dtc, D_table, dz, K_table, m, mob_fak, position_z, ths, v, z):
    # Initialisation of parameters and arrays

    position_znew = copy.deepcopy(position_z) # position of every particle in soil matrix
    particles = np.zeros((dim,1)) # amount of particles in each grid element of soil matrix
    theta = np.zeros((dim,1)) # soil moisture profile of soil matrix
    Cw = np.zeros((dim,1)) # concentration profile of soil matrix
    mtx_avg_age = np.zeros((dim,1)) # age profile of soil matrix  
    M_layer = np.zeros((dim,1)) # total water mass in each layer
    ip = np.zeros((dim,1)) # storage for particle amounts in each grid element

    ip_0 = position_z[:,0] > z[0] # correction for unphysical particle losses at top of the soil matrix
    position_z[ip_0,0] = z[0] - position_z[ip_0,0] # ensures that every particle has a negative position value within the soil matrix

    # Particle displacement
    for j in range(0,dim-1):
        icl = np.sum(D_table < D[j])  # finds the last bin where the assumed diffusivity of the lookup table (D_table) is still smaller than the real diffusivity (D). 'icl' is then the number of the bin with the maximum possible diffusivity 
        ip[j+1] = np.sum(np.logical_and(position_z[:,0] <= z[j],position_z[:,0] > z[j+1])) # finds the amount of all particles within the actual grid element (between z(j) and z(j+1))
        ip_count = np.array(np.sum(np.logical_and(position_z[:,0] <= z[j],position_z[:,0] > z[j+1])))
    
        if ip[j+1] > 0:
            d_part = np.floor(ip[j+1] / (icl+1)) # allocates the total amount of particles within a grid element to the possible amount of diffussivity bins. 'd_part' is then the total amount of particles within a diffussivity bin
            mobile = np.floor(mob_fak * ip[j+1]) # fraction of particles which contribute to the particle movement, thus the "mobile" fraction
     
            # defining of scaling factors
            mo_fak = np.zeros((ip_count,1)) # initialisation of factor for scaling diffussive motion of particles in different bins
            k_fak = np.zeros((ip_count,1)) # initialisation of factor for scaling advective motion of particles in different bins
        
            for ii in range(0,icl): # establishes a factor for scaling the diffusivity of the particles stored in the different diffussivity bins
                mo_fak[(ii * (d_part[0,].astype(int))+1):(ii+1) * d_part[0,].astype(int)+1] = (D_table[0,ii] / D[j]) * np.ones((d_part[0,].astype(int),1)) # 'mo_fak' is thereby just the relation of the diffusivity of the single bins to the real diffussivity in the actual grid element. This factor is then applied to every particle within the bins
                k_fak[(ii * (d_part[0,].astype(int)))+1:(ii+1) * d_part[0,].astype(int)+1] = (K_table[0,ii] / v[j]) * np.ones((d_part[0,].astype(int),1)) # 'k_fak' is thereby just the relation of the hydraulic conductivity of the single bins to the real hydraulic conductivity in the actual grid element. This factor is then applied to every particle within the bins
     
            mo_fak[(icl * (d_part[0,].astype(int))):(mo_fak.size)] = np.ones((mo_fak.size-(icl*d_part[0,].astype(int)),1)) # these particles not considered by the floor function just get the highest scaling factor of one
            k_fak[(icl * (d_part[0,].astype(int))):(k_fak.size)] = np.ones((k_fak.size-((icl)*d_part[0,].astype(int)),1)) # these particles not considered by the floor function just get the highest scaling factor of one

            # interpolates v and D for all particles within the grid element
            zint = (np.absolute(position_z[(np.sum(ip)-ip[j+1]).astype(int),0]) - np.absolute(z[j])) / (np.absolute(z[j+1]) - np.absolute(z[j]))
            v_part = v[j] - (zint * (v[j]-v[j+1]))
            D_part = D[j] - (zint * (D[j]-D[j+1]))

            # calculates the random walk term
            correction = 0.25 * (D[j]-D[j+1]) * mo_fak[(ip_count-mobile[0,].astype(int)):(ip_count),0] / dz[0] # correction term for random walk
            randomwalk = 1 * 1 * (np.random.randn(mobile[0,].astype(int),1))[:,0] * mo_fak[(ip_count-mobile[0,].astype(int)):ip_count,0] *  np.sqrt(2 * D_part * dtc)  # random walk term
            position_znew[(np.sum(ip).astype(int)-mobile[0,].astype(int)):np.sum(ip).astype(int),0] = position_z[(np.sum(ip).astype(int)-mobile[0,].astype(int)):np.sum(ip).astype(int),0] - ((1 * v_part * k_fak[ip_count-mobile[0,].astype(int):ip_count,0]+correction) * dtc + randomwalk[:,]) # calculates the new positions of the particles in each grid element
                       
    # postprocessing for realistic particle positions
    ip_0 = position_znew[:,0] >= z[0] # correction for unphysical particle losses at top of the soil matrix due to the displacement process
    position_znew[ip_0,0] = z[0] - position_znew[ip_0,0]

    ip_low = position_znew[:,0] <= z[dim-1] # correction for unphysical particle losses at bottom of the soil matrix due to the displacement process
    position_znew[ip_low,0] = z[dim-2] 


    #Calculation of water and solute masses in each grid element and generating of new soil moisture and concentration profile
    for i in range (0,dim-1):
        ipart = np.logical_and(position_znew[:,0] <= z[i],position_znew[:,0] > z[i+1])
        particles[i] = np.sum(ipart) # amount of particles in grid element
    
        # check for oversaturation
        if particles[i] > np.round(ths * dz[i] * 1 * 1000 / m): # if the number of particles is higher than the saturated threshold the number is set to this threshold
            particles[i] = np.round(ths * dz[i] * 1 * 1000 / m)
       
            ipart2 = np.where(np.logical_and(position_znew[:,0] <= z[i],position_znew[:,0] > z[i+1]))
            for ii in range(np.round(ths * dz[i] * 1 * 1000 / m)+1,ipart2[0][:].size): 
                 position_znew[ipart2[0][ii]]= position_znew[ipart2[0][ii]] - dz[i]
    
        # calculation of new profiles
        M_layer[i] = particles[i] * m # total mass in each grid element
        theta[i] = np.round(M_layer[i] / (dz[i]*1*1000),4) # new soil moisture in each grid element
    
        mtx_avg_age[i] = np.sum(position_znew[ipart,1]) / particles[i] # new average age in each grid element
        if np.isnan(mtx_avg_age[i]):
            mtx_avg_age[i] = 0
        
        Cw[i]= np.sum(position_znew[ipart,2]) / particles[i] # new concentration in each grid element
        position_znew[ipart,2] = Cw[i] * np.ones((np.sum(ipart),)) # new concentration of every single particle in each grid element
    
    return Cw, mtx_avg_age, position_znew, theta

## displacement routine for event particles in matrix
def displ_mtx_event(dim, dtc, D_table, K_table, m, dz, position_z_event, prob, ths, z):
    # Initialisation of parameters and arrays

    position_znew_event = copy.deepcopy(position_z_event) # position of every event particle in soil matrix
    particles = np.zeros((dim,1)) # amount of particles in each grid element of soil matrix
    theta_event = np.zeros((dim,1)) # event soil moisture profile of soil matrix in the first grid elements
    M_layer = np.zeros((dim,1)) # total water mass in each layer
    ip = np.zeros((dim,1)) # storage for particle amounts in each grid element

    ip_0 = position_z_event[:,0] > z[0] # correction for unphysical particle losses at top of the soil matrix
    position_z_event[ip_0,0] = z[0] - position_z_event[ip_0,0] # ensures that every particle has a negative position value within the soil matrix

    # Particle displacement
    for j in range(0,dim-1): 
       ip[j+1] = np.sum(np.logical_and(position_z_event[:,0] <= z[j],position_z_event[:,0] > z[j+1])) # finds the amount of all particles within the actual grid element (between z(j) and z(j+1))
       ip_count = np.array(np.sum(np.logical_and(position_z_event[:,0] <= z[j],position_z_event[:,0] > z[j+1])))
       
       if ip[j+1] > 0: 
          # calculation of v and D
          v_part = 2.5 * K_table[0,K_table.size-1] # takes the highest K value for calculating advective velocity
          D_part = np.quantile(D_table,prob) # takes the highest D value depending on the "prob-quantile" for calculating diffussivity
          
          # calculation of random walk term
          correction = 0 # correction term for random walk 
          randomwalk = 1 * 1 * (2 * np.random.rand(ip_count,1) - np.ones((ip_count,1))) * np.sqrt(6 * D_part * dtc) # random walk term ; verteilung zwischen -1 und 1 (gleichverteilt)
          position_znew_event[ip[j,0].astype(int):np.sum(ip).astype(int),0] = position_z_event[ip[j,0].astype(int):np.sum(ip).astype(int),0] - ((v_part + correction) * np.ones((ip_count,1))[:,0] * dtc + randomwalk[:,0]) # calculates the new positions of event water particles for the first grid elements
       
      
    # postprocessing for realistic particle positions
    ip_0 = position_znew_event[:,0] > z[0] # correction for unphysical particle losses at top of the soil matrix due to the displacement process
    position_znew_event[ip_0,0] = z[0] - position_znew_event[ip_0,0]

    ip_low = position_znew_event[:,0] <= z[dim-1] # correction for unphysical particle losses at bottom of the soil matrix due to the displacement process
    position_znew_event[ip_low,0] = z[dim-2] 

    # Calculation of water and solute masses in each grid element and generating of new soil moisture and concentration profile

    for i in range(0,dim-1):
        ipart = np.logical_and(position_znew_event[:,0] <= z[i],position_znew_event[:,0] > z[i+1])
        particles[i] = np.sum(ipart) # amount of particles in grid element
    
        # check oversaturation
        if particles[i] > np.round(ths *dz[i] * 1 * 1000 / m): # if the amount of particles is higher than the saturated threshold, the amount is set to this threshold
            particles[i] = np.round(ths *dz[i] * 1 * 1000 / m)
       
            ipart2 = np.where(np.logical_and(position_znew_event[:,0] <= z[i],position_znew_event[:,0] > z[i+1]))
            for ii in range(np.round(ths * dz[i] * 1 * 1000 / m)+1,ipart2[0][:].size): 
                 position_znew_event[ipart2[0][ii]]= position_znew_event[ipart2[0][ii]] - dz[i]
        
        # calculation of new event moisture profile 
        M_layer[i] = particles[i] * m # total mass in each layer (kg)
        theta_event[i] = M_layer[i] / (dz[i] * 1 * 1000) # soil moisture in each layer (kg/m^3)
   
        val_pos = -np.sort(-position_znew_event[:,0])
        idx = np.argsort(-position_znew_event[:,0])
        
        position_znew_event = np.concatenate((np.reshape(val_pos, (val_pos[:,].size,1)),np.reshape(position_znew_event[idx,1],(position_znew_event[idx,1].size,1)),np.reshape(position_znew_event[idx,2],(position_znew_event[idx,2].size,1)),np.reshape(position_znew_event[idx,3],(position_znew_event[idx,3].size,1))), axis=1)               

    return position_znew_event, theta_event

## displacement routine for particles in pfd
def displ_pfd(n_mak, pfd_dim, pfd_dz, pfd_m, pfd_n, pfd_position_z, pfd_r, pfd_z):
    
    # Initialisation of parameters and arrays
    pfd_position_znew = copy.deepcopy(pfd_position_z) # position of every particle in pfd
    pfd_particles = np.zeros((pfd_dim,1)) # amount of pfd_particles in each grid element of pfd
    pfd_theta = np.zeros((pfd_dim,1)) # soil moisture profile of pfd
    pfd_Cw = np.zeros((pfd_dim,1)) #concentration profile of pfd
    avg_age = np.zeros((pfd_dim,1)) # average age progile of pfd (only for testin, yet)
    M_layer = np.zeros((pfd_dim,1)) # total mass in each layer
    ip = np.zeros((pfd_dim,1))

    ip_0 = pfd_position_z[:,0] > pfd_z[0] # correction for unphysical particle losses at top of the pfd
    pfd_position_z[ip_0,0] = pfd_z[0] - pfd_position_z[ip_0,0] # ensures that every particle has a negative position value within the pfd

    # Particle displacement
    for j in range(0,pfd_dim-1): 
        ip[j+1] = np.sum(np.logical_and(pfd_position_z[:,0] <= pfd_z[j],pfd_position_z[:,0] > pfd_z[j+1])) # finds all pfd_particles within the actual grid element (between z(j) and z(j+1))
       # ip_count = np.array(np.sum(np.logical_and(pfd_position_z[:,0] <= pfd_z[j],pfd_position_z[:,0] > pfd_z[j+1])))
    
        if ip.size != 0:
           pfd_position_znew[ip[j,0].astype(int):np.sum(ip).astype(int),0] = pfd_z[pfd_dim-1] # every particle in the respective grid elements is initially displaced to the last grid element/bottom of the macropore
            

    # Calculation of water and solute masses in each grid element and generating of new soil moisture and concentration profile
    for i in range(pfd_dim-1,0,-1):
        ipart = np.where(np.logical_and(pfd_position_znew[:,0] < pfd_z[i-1],pfd_position_znew[:,0] >= pfd_z[i]))
        pfd_particles[i] = ipart[:][0].size # amount of pfd_particles in grid element
    
        # check oversaturation
        if pfd_particles[i] > np.round(pfd_n[i-1] * n_mak): # if the number of pfd_particles is higher than the saturated threshold the number is set to this threshold
            pfd_particles[i] = np.round(pfd_n[i-1] * n_mak)
       
            for ii in range(ipart[0][0].astype(int),(ipart[0][-1] - np.round(pfd_n[i-1] * n_mak)).astype(int)): # redistribution of residual pfd_particles out of the last grid element to the residual unsaturated grid elements above
                pfd_position_znew[ii,0] = pfd_position_znew[ii,0] + pfd_dz[i-1] 
        
        ip_0 = pfd_position_znew[:,0] > pfd_z[0] # again correction for unphysical particle losses at top of the pfd
        pfd_position_znew[ip_0,0] = pfd_z[0]  
   
        blabla = np.where(pfd_position_znew[:,0] == np.round(pfd_z[i],2))
        pfd_position_znew[blabla[:][0],0] = pfd_z[i] + (np.round(pfd_z[i-1] - pfd_z[i],2) * np.random.rand(blabla[:][0].size,)) # every particle gets randomly an own position in the respective grid element (just for style)
    
        # calculation of new profiles
        M_layer[i] = pfd_particles[i] * pfd_m # total mass in each grid element
        pfd_theta[i] = np.round(M_layer[i] / ((pfd_dz[i-1] * (np.pi * pfd_r**2) * 1000) *n_mak),3) # new soil moisture in each grid element

        ipart2 = np.where(np.logical_and(pfd_position_znew[:,0] < pfd_z[i-1],pfd_position_znew[:,0] >= pfd_z[i]))
        avg_age[i] = np.sum(pfd_position_znew[ipart2[:][0],1]) / pfd_particles[i] # new average age in each grid element
        pfd_Cw[i] = np.sum(pfd_position_znew[ipart2[:][0],2]) / pfd_particles[i] # new concentration in each grid element
    
        if np.isnan(pfd_Cw[i]):
            pfd_Cw[i] = 0
        
        if np.isnan(avg_age[i]):
            avg_age[i] = 0
        
    
        pfd_position_znew[ipart2[:][0],2] = pfd_Cw[i] * np.ones((ipart2[:][0].size,)) # new concentration of every single particle in each grid element
    
    return pfd_particles, pfd_Cw, pfd_theta, pfd_position_znew

## mixing routine for mixing of pre-event particles with event particles in the matrix
def mixing_mtx(age_event, Cw_event, position_z, position_z_event, time):
    
    ip_mix = np.where(position_z_event[:,1] < (time - position_z_event[:,3])) # finds the positions of all particles which are old enough to contribute to the mixing process at given mixing time and elapsed simulation time
    ip_res = position_z_event[:,1] >= (time - position_z_event[:,3]) # finds the positions of the residual particles which are still to young for mixing
   
    if ip_mix[:][0].size > 0:
        ipp = np.where(ip_mix[:][0] < position_z_event[:,0].size) # selects from all particles just these which are able to contribute to mixing
        
        #adds the mixing particles to the position and concentration arrays of the pre-event particles and deletes them out of the event arrays                
#        pos_age_conc_mix = pd.concat([pd.DataFrame(position_z_event[ip_mix[0][ipp[:][0].astype(int)],0]), pd.DataFrame(position_z_event[ip_mix[0][ipp[:][0].astype(int)],1]), pd.DataFrame(position_z_event[ip_mix[0][ipp[:][0].astype(int)],2])], axis=1)
#        pos_age_conc_mix = np.asarray(pos_age_conc_mix)
        
        pos = np.reshape(position_z_event[ip_mix[0][ipp[:][0]],0], (-1,1))
        age = np.reshape(position_z_event[ip_mix[0][ipp[:][0]],1], (-1,1))
        conc = np.reshape(position_z_event[ip_mix[0][ipp[:][0]],2], (-1,1))
        
        pos_age_conc_mix = np.concatenate((pos, age, conc), axis = 1)
        
        position_z = np.append(pos_age_conc_mix,position_z, axis=0)
                          
        position_z_event = position_z_event[ip_res,:]
        
    return age_event, Cw_event, position_z, position_z_event

## mixing routine for mixing of particles from pfd into matrix
def mixing_pfd_mtx(dtc, k, ks, m, mak_sml, mak_mid, n_mak, particle_contact_grid, pfd_dim, pfd_dz, pfd_m, pfd_n, pfd_particles, pfd_position_z, pfd_r, pfd_theta, pfd_z, position_z, psi, rate_big, rate_mid, rate_sml, rho, theta, ths, z):
    
    for yy in range(pfd_dim-1,0,-1):  

        # Finds the respective depths of the soil matrix corresponding to the depths of the three macropore sizes
        indx_big = np.where(pfd_z[yy-1] > np.round(z,2)) # for big macropore

        if ((yy-1)-((pfd_dim-1) - np.where(np.round(pfd_z,2) == mak_mid)[:][0])-1) > 0: # for medium macropores
           indx_mid = np.where(pfd_z[(yy-1) - ((pfd_dim-1) - np.where(np.round(pfd_z,2) == mak_mid)[:][0])-1] > np.round(z,2))
        else:
           indx_mid = np.nan
            

        if ((yy-1)-((pfd_dim-1) - np.where(np.round(pfd_z,2) == mak_sml)[:][0])-1) > 0: # for small macropores
            indx_sml = np.where(pfd_z[(yy-1) - ((pfd_dim-1) - np.where(np.round(pfd_z,2) == mak_sml)[:][0])-1] > np.round(z,2))
        else:
            indx_sml = np.nan
            
        # Calculates mixing mass and converts it into integer number of particles entering the soil matrix
        # assumption: diffusive mixinng only from saturated grid elements of pfd into soil matrix, no back diffusion from matrix into pfd
        if np.logical_and(pfd_theta[yy] == 1,(theta[indx_big[0][0]] + theta[indx_big[0][0]-1])/2 < ths):

            # calculates the flux density and the mass of diffusive mixing based on psi and converts the mass to integer number of particles leaving the pfd 
            q_mix_psi = np.abs((2 * ks * ((k[indx_big[0][0]] + k[indx_big[0][0]-1]) / 2)) / (ks + ((k[indx_big[0][0]] + k[indx_big[0][0]-1]) / 2)) * (((psi[indx_big[0][0]] + psi[indx_big[0][0]-1]) / 2) / (pfd_r*2))) # diffussive flow and mass calculated with psi_init
            m_mix_psi = (q_mix_psi * (2 * np.pi * pfd_r * np.round(pfd_dz[yy-1],2))) * rho * dtc * n_mak # calculation of total mass leaving pfd in current depth
            pfd_particles_mix = np.floor(m_mix_psi / pfd_m) # amount particles leaving macropore

            # calculates the number of particles having contact to the lateral surface of a pfd grid element and ensures that only the maximum possible number of contact particles is leaving
            rate = (particle_contact_grid * n_mak) / (pfd_n[yy-1] * n_mak) 
            pfd_particles_contact = np.floor(pfd_particles[yy] * rate) 

            if pfd_particles_mix > pfd_particles_contact: # maximum possible contact particles
               pfd_particles_mix = copy.deepcopy(pfd_particles_contact)

            # converts leaving pfd particles into integer number of matrix particles entering the soil matrix
            mtx_particles_mix = np.floor((pfd_particles_mix * pfd_m) / m)

            pfd_particles_mix = np.round(pfd_particles_mix - ((((pfd_particles_mix * pfd_m) / m) - mtx_particles_mix) * m) / pfd_m) # residual amount of particles not considered by "floor" function are recalculated to the pfd            
           
            indx2 = np.where(np.logical_and(pfd_position_z[:,0] < pfd_z[yy-1], pfd_position_z[:,0] >= pfd_z[yy]))            
            pfd_age_output = np.sum(pfd_position_z[(indx2[0][-1]-pfd_particles_mix[0,].astype(int)):indx2[0][-1],1]) / pfd_particles_mix[0,]
            pfd_conc_output2 = np.sum(pfd_position_z[(indx2[0][-1]-pfd_particles_mix[0,].astype(int)):indx2[0][-1],2]) / pfd_particles_mix[0,] # solute concentration of leaving particles
            
            # Deletion of leaving particle amount out of the pfd arrays
            pfd_position_z = np.delete(pfd_position_z, np.array(range((indx2[0][-1]-pfd_particles_mix[0,].astype(int)),indx2[0][-1])) , axis =0)
            
            # Adds the amount of particle entering the soil matrix to the arrays in the respective depths
            # determines how the total number of entering particles is distributed to  the different soil matrix depths based on the actual grid elements of the three macropore sizes
            
            if np.isnan(indx_sml).any() == False: # valid when all macropores are still unsaturated
               pos_sml = np.round(z[indx_sml[0][0]] + np.round(z[indx_sml[0][0]-1] - z[indx_sml[0][0]],1) * np.random.rand((np.round(mtx_particles_mix*rate_sml))[0,].astype(int),1),3)
               pos_mid = np.round(z[indx_mid[0][0]] + np.round(z[indx_mid[0][0]-1] - z[indx_mid[0][0]],1) * np.random.rand((np.round(mtx_particles_mix*rate_mid))[0,].astype(int),1),3)
               pos_big = np.round(z[indx_big[0][0]] + np.round(z[indx_big[0][0]-1] - z[indx_big[0][0]],1) * np.random.rand((np.round(mtx_particles_mix*rate_big))[0,].astype(int),1),3)
               
               pos = np.concatenate((pos_big, pos_mid, pos_sml))
               age = np.concatenate((pfd_age_output * np.ones((pos_big.size,1))[:,], pfd_age_output * np.ones((pos_mid.size,1))[:,], pfd_age_output * np.ones((pos_sml.size,1))[:,]))
               conc = np.concatenate((pfd_conc_output2 * np.ones((pos_big.size,1))[:,], pfd_conc_output2 * np.ones((pos_mid.size,1))[:,], pfd_conc_output2 * np.ones((pos_sml.size,1))[:,]))
               pos_age_conc_mix = np.concatenate((pos, age, conc), axis=1)
               
               position_z = np.concatenate((pos_age_conc_mix,position_z), axis=0)
               
               val_pos = -np.sort(-position_z[:,0])
               idx = np.argsort(-position_z[:,0])
               
               position_z = np.concatenate((np.reshape(val_pos, (val_pos[:,].size,1)),np.reshape(position_z[idx,1],(position_z[idx,1].size,1)),np.reshape(position_z[idx,2],(position_z[idx,2].size,1))), axis=1)
           
            elif np.isnan(indx_mid).any() == False: # valid when just big and medium macropores are still unsaturated
               pos_big = np.round(z[indx_big[0][0]] + np.round(z[indx_big[0][0]-1] - z[indx_big[0][0]],1) * np.random.rand((np.floor(mtx_particles_mix * ((rate_sml / 2) + rate_big)))[0,].astype(int),1),3)
               pos_mid = np.round(z[indx_mid[0][0]] + np.round(z[indx_mid[0][0]-1] - z[indx_mid[0][0]],1) * np.random.rand((np.ceil(mtx_particles_mix * ((rate_sml / 2) + rate_mid)))[0,].astype(int),1),3)
            
               pos = np.concatenate((pos_big, pos_mid))
               age = np.concatenate((pfd_age_output * np.ones((pos_big.size,1))[:,], pfd_age_output * np.ones((pos_mid.size,1))[:,]))
               conc = np.concatenate((pfd_conc_output2 * np.ones((pos_big.size,1))[:,], pfd_conc_output2 * np.ones((pos_mid.size,1))[:,]))
               pos_age_conc_mix = np.concatenate((pos, age, conc), axis=1)
                
               position_z = np.concatenate((pos_age_conc_mix,position_z), axis=0)
               
               val_pos = -np.sort(-position_z[:,0])
               idx = np.argsort(-position_z[:,0])
               
               position_z = np.concatenate((np.reshape(val_pos, (val_pos[:,].size,1)),np.reshape(position_z[idx,1],(position_z[idx,1].size,1)),np.reshape(position_z[idx,2],(position_z[idx,2].size,1))), axis=1)               
               
            else: # valid when just big macropores are still unsaturated
               pos_big = np.round(z[indx_big[0][0]] + np.round(z[indx_big[0][0]-1] - z[indx_big[0][0]],1) * np.random.rand((np.round(mtx_particles_mix * (rate_sml + rate_mid + rate_big)))[0,].astype(int),1),3)
            
               pos = pos_big
               age = pfd_age_output * np.ones((pos_big.size,1))[:,]
               conc = pfd_conc_output2 * np.ones((pos_big.size,1))[:,]
               pos_age_conc_mix = np.concatenate((pos, age, conc), axis=1)
               
               position_z = np.concatenate((pos_age_conc_mix,position_z), axis=0)
               
               val_pos = -np.sort(-position_z[:,0])
               idx = np.argsort(-position_z[:,0])
               
               position_z = np.concatenate((np.reshape(val_pos, (val_pos[:,].size,1)),np.reshape(position_z[idx,1],(position_z[idx,1].size,1)),np.reshape(position_z[idx,2],(position_z[idx,2].size,1))), axis=1)                                             
    
    return pfd_position_z, position_z


        
        
  
            
    
            