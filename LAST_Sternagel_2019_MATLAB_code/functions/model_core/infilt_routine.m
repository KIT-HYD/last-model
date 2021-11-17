function [age_event, bla, bla1, Cw_event, mass_totalsoluteinput, mass_totalwaterinput, m_input, m_surface, m_slt_surface, pfd_age, pfd_Cw_event,...
         pfd_position_z, position_z_event, solute_test, solute_test2, theta] =... % output arguments
            infilt_routine... % function
            (age_event, bla, bla1, Cw_event, Cw_eventparticles, dtc, dz, k, ks, m, mass_totalsoluteinput, mass_totalwaterinput, m_input,... % input arguments
            m_slt_surface, m_surface, n_mak, pfd_age, pfd_Cw_event, pfd_dz, pfd_k, pfd_m, pfd_position_z, pfd_r, pfd_theta, pfd_z,...
            position_z_event, psi_init, qb_u, rho, solute_test, solute_test2, theta, ths, time, z)

             
        
%% Infiltration routine defining different scenarios 
%% for distribution of the precipitation water to the soil matrix and the pfd

if (abs(qb_u)> 0) || (m_surface > 0)

%% Initialisation of input masses and surface storages

 m_surface = m_input + ((abs(qb_u)*1000)*dtc) + m_surface; % rainwater initially stored in a virtual water surface storage
 m_input = 0;
 m_slt_surface = m_slt_surface + (((abs(qb_u)*1000)*dtc)*Cw_eventparticles); % solute mass initially stored in a virtual solute surface storage
 
 m_input1 = ((k(1)+ks)/2) * (((abs(psi_init(1)-psi_init(2)))/dz(1))+1) * 1 * rho * dtc; % water mass input into soil matrix
 m_input2 = ((pfd_k*((pi*pfd_r^2)))*rho*dtc) * n_mak ; % water mass input into pfd

 mass_totalwaterinput = mass_totalwaterinput + (abs(qb_u)*1000) * dtc; % control parameter, sums up the total incoming water mass
 mass_totalsoluteinput = mass_totalsoluteinput + (((abs(qb_u)*1000)*dtc) * Cw_eventparticles); % control parameter, sums up the total incoming solute mass

     if Cw_eventparticles == 0 % concentration input after rainfall event dependent on masses in surface storages
        Cw_eventparticles = m_slt_surface/m_surface; 
     end

 theta_count = 0; % scenario counter    
     
 %% Scenario 1: Initially water and solutes flow into the soil matrix and residual overflow infiltrates simultaneously the pfd
 
 if theta(1) < (ths*0.8)
     
     if m_input1 > m_surface
        m_input1 = m_surface;
     end
     
     if Cw_eventparticles*m_input1 > m_slt_surface
        Cw_eventparticles = m_slt_surface/m_input1; 
     end
     
     % calls function describing infiltration into soil matrix
     [age_event, bla, Cw_event, m_input, m_slt_surface, m_surface, position_z_event, solute_test, theta, theta_count] = ... % output arguments
        infilt_mtx... % function
            (age_event, bla, Cw_event, Cw_eventparticles, dz, m, m_input1, m_slt_surface, m_surface, position_z_event, ... % input arguments
            rho, solute_test, theta, time, z); 
     
 end

 % if the incoming precipitation mass is higher than the infiltrating mass
 % (m_input1) it comes to a water accumulation of this residual water at the surface.
 % This accumulation then flows immediately into the pfd with m_input2.
 
 if m_surface > 0
     
     if m_input2 > m_surface
        m_input2 = m_surface;
     end

     % calls function describing infiltration into pfd
     [bla1, m_input, m_slt_surface, m_surface, pfd_age, pfd_Cw_event, pfd_position_z, pfd_theta, solute_test2, theta_count] = ... % output arguments
        infilt_pfd... % function
            (bla1, Cw_eventparticles, m_input2, m_slt_surface, m_surface, n_mak, pfd_age, pfd_Cw_event, pfd_dz, ... % input arguments
            pfd_m, pfd_position_z, pfd_r, pfd_theta, pfd_z, rho, solute_test2, theta_count,  time);
    
 end
 
%% Scenario 2: If the top grid element of the soil matris is saturated the water and solutes infiltrate only the pfd

if theta_count == 0 && round(theta(1),2) == ths && (pfd_theta(1) < 1)
    
    if m_input2 > m_surface
       m_input2 = m_surface;
    end
    
    % calls function describing infiltration into pfd
    [bla1, m_input, m_slt_surface, m_surface, pfd_age, pfd_Cw_event, pfd_position_z, pfd_theta, solute_test2, theta_count] = ... % output arguments
        infilt_pfd... % function
            (bla1, Cw_eventparticles, m_input2, m_slt_surface, m_surface, n_mak, pfd_age, pfd_Cw_event, pfd_dz, ... % input arguments
            pfd_m, pfd_position_z, pfd_r, pfd_theta, pfd_z, rho, solute_test2, theta_count,  time);
        
end 


end


end

