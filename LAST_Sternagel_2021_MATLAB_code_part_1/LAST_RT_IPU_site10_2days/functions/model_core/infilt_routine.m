function [bla, bla1, mass_totalsoluteinput, mass_totalwaterinput, m_surface, m_slt_surface,...
         pfd_position_z, position_z_event, solute_test, solute_test2, theta] =... % output arguments
            infilt_routine... % function
            (istyp, reactive_transport, pfd_dim,t_mix_mean2,t_mix_SD2,t_mix_mean1,t_mix_SD1,t_mix_ratio, bla, bla1, Cw_eventparticles, dtc, dz, k, ks, m, mass_totalsoluteinput, mass_totalwaterinput,... % input arguments
            m_slt_surface, m_surface, n_mak, pfd_dz, pfd_k, pfd_m, pfd_position_z, pfd_r, pfd_theta, pfd_z,...
            position_z_event, psi, qb_u, rho, solute_test, solute_test2, theta, ths, time, z)

             
        
%% Infiltration routine defining different scenarios 
%% for distribution of the precipitation water to the soil matrix and the pfd

if (abs(qb_u)> 0) || (m_surface > 0)

%% Initialisation of input masses and surface storages

 m_surface = ((abs(qb_u)*1000*1)*dtc) + m_surface; % rainwater initially stored in a virtual water surface storage
 m_slt_surface = m_slt_surface + (((abs(qb_u))*1*dtc)*Cw_eventparticles); % solute mass initially stored in a virtual solute surface storage
 Cw_eventparticles = m_slt_surface / (m_surface/rho);
 
 m_input1 = ((k(1)+ks(istyp(1)))/2) * (((abs(0-psi(1)))/dz(1))+1) * 1 * rho * dtc; % water mass input into soil matrix
 m_input2 = ((pfd_k*((pi*pfd_r^2)))*rho*dtc) * n_mak ; % water mass input into pfd
 

 mass_totalwaterinput = mass_totalwaterinput + (abs(qb_u)*1000) * dtc; % control parameter, sums up the total incoming water mass
 mass_totalsoluteinput = mass_totalsoluteinput + (((abs(qb_u))*dtc) * Cw_eventparticles); % control parameter, sums up the total incoming solute mass

     if Cw_eventparticles == 0 % concentration input after rainfall event dependent on masses in surface storages
        Cw_eventparticles = m_slt_surface/(m_surface/rho); 
     end

 theta_count = 0; % scenario counter    
     
 %% Scenario 1: Initially water and solutes flow into the soil matrix and residual overflow infiltrates simultaneously the pfd
 
 % if the incoming precipitation mass is higher than the infiltrating mass
 % (m_input1) it comes to a water accumulation of this residual water at the surface.
 % This accumulation then flows immediately into the pfd with m_input2.
 
 if m_surface > 0
     
     if m_input2 > m_surface
        m_input2 = m_surface;
     end

     % calls function describing infiltration into pfd
     [bla1, m_slt_surface, m_surface, pfd_position_z, pfd_theta, solute_test2, theta_count] = ... % output arguments
        infilt_pfd... % function
            (reactive_transport, pfd_dim,t_mix_mean2,t_mix_SD2,t_mix_mean1,t_mix_SD1,t_mix_ratio, bla1, Cw_eventparticles, m_input2, m_slt_surface, m_surface, n_mak, pfd_dz, ... % input arguments
            pfd_m, pfd_position_z, pfd_r, pfd_theta, pfd_z, rho, solute_test2, theta_count,  time);
    
 end

 if theta(1) < (ths(1) * 1.0)
     
     if m_input1 > m_surface
        m_input1 = m_surface;
     end
     
     if Cw_eventparticles*(m_input1/rho) > m_slt_surface
        Cw_eventparticles = m_slt_surface/(m_input1/rho); 
     end
     
     % calls function describing infiltration into soil matrix
     if m_input1 > m
     [bla, m_slt_surface, m_surface, position_z_event, solute_test, theta, theta_count] = ... % output arguments
        infilt_mtx... % function
            (t_mix_mean2,t_mix_SD2,t_mix_mean1,t_mix_SD1,t_mix_ratio, bla, Cw_eventparticles, dz, m, m_input1, m_slt_surface, m_surface, position_z_event, ... % input arguments
            rho, solute_test, theta, time, z); 
     end
 end

 
%% Scenario 2: If the top grid element of the soil matris is saturated the water and solutes infiltrate only the pfd

if theta_count == 0 && round(theta(1),2) == ths(1) && (pfd_theta(1) < 1)
    
    if m_input2 > m_surface
       m_input2 = m_surface;
    end
    
    % calls function describing infiltration into pfd
     [bla1, m_slt_surface, m_surface, pfd_position_z, pfd_theta, solute_test2, theta_count] = ... % output arguments
        infilt_pfd... % function
            (reactive_transport, pfd_dim,t_mix_mean2,t_mix_SD2,t_mix_mean1,t_mix_SD1,t_mix_ratio, bla1, Cw_eventparticles, m_input2, m_slt_surface, m_surface, n_mak, pfd_dz, ... % input arguments
            pfd_m, pfd_position_z, pfd_r, pfd_theta, pfd_z, rho, solute_test2, theta_count,  time);
end 


end


end

