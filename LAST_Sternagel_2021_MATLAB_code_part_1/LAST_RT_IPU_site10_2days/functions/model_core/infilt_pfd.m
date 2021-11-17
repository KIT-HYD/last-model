%% Function to describe the infiltration into pfd

function [bla1, m_slt_surface, m_surface, pfd_position_z, pfd_theta, solute_test2, theta_count] = ... % output arguments
            infilt_pfd... % function
                (reactive_transport, pfd_dim,t_mix_mean2,t_mix_SD2,t_mix_mean1,t_mix_SD1,t_mix_ratio, bla1, Cw_eventparticles, m_input2, m_slt_surface, m_surface, n_mak, pfd_dz, ... % input arguments
                pfd_m, pfd_position_z, pfd_r, pfd_theta, pfd_z, rho, solute_test2, theta_count,  time)


for i = 1:pfd_dim-1
if (m_input2 > pfd_m) && (pfd_theta(i) < 1)

m_input3 = m_input2;   
pfd_theta_diff = 1 - pfd_theta(i); % detects the amount of water the first grid element of the pfd can gather in this time step
precip_particles = floor((m_input3)/pfd_m); % amount of event particles to be injected, converts the mass input to number of particles
m_surface = round((m_surface - m_input3) + (m_input3-(precip_particles*pfd_m)),4) ; % updates water surface storage
theta_count = 1; % sets scenario counter to 1, thus the routine will not run through the subsequent scenarios
pfd_theta(i)=pfd_theta(i)+ (precip_particles*pfd_m)/((rho*pfd_dz(1)*pi*(pfd_r^2))*n_mak); % updates the soil moisture of the first grid element

        % if the first grid element of the pfd is oversaturated due to a
        % too high input mass (m_input2) then the surplus is considered by
        % updating the following parameters
        if pfd_theta(i) > 1
           m_surface = m_surface + (pfd_theta(i)-1)*((rho*pfd_dz(1)*pi*(pfd_r^2))*n_mak);
           precip_particles = floor((pfd_theta_diff*((rho*pfd_dz(1)*pi*(pfd_r^2))*n_mak)/pfd_m));
           m_input3 = precip_particles * pfd_m;
           pfd_theta(i) = 1;
        end
        
        
        % ensures that the infiltrating solute mass is not greater than the
        % solute mass storage
        if Cw_eventparticles*(m_input3/rho) > m_slt_surface
           Cw_eventparticles = m_slt_surface/(m_input3/rho); 
        end

m_input2 = m_input2 - m_input3;
m_slt_surface = (m_slt_surface - (Cw_eventparticles*(m_input3/rho))) + (Cw_eventparticles*((m_input3-(precip_particles*pfd_m))/rho)); % updates solute mass surface storage
solute_test2 = solute_test2 + (Cw_eventparticles*((precip_particles*pfd_m)/rho)); % counter for amount of solutes totally infiltrating the soil matrix within the first scenario

bla1 = bla1 + (precip_particles*pfd_m); % counter for amount of water totally infiltrating the pfd
pfd_precip_position = pfd_z(i) * ones(precip_particles,1); % in every time step the new infiltrating particles get the position pfd_z(1)
pfd_age = time * ones(precip_particles,1); % matrix with the age of the infiltrating event particles, every new particle gets the value of its timestep when entering the domain

%adds the number of event particles entering the pfd with
%their concentration value to the pfd_Cw_event matrix
if Cw_eventparticles == 0                    
   pfd_Cw_event = zeros(precip_particles,1);
else
   pfd_Cw_event = Cw_eventparticles*ones(precip_particles,1);
end

pfd_t_mix = [randn(floor(precip_particles*t_mix_ratio),1) * t_mix_SD1 + t_mix_mean1;randn(ceil(precip_particles*(1-t_mix_ratio)),1) * t_mix_SD2 + t_mix_mean2];
pfd_t_mix = pfd_t_mix(randperm(numel(pfd_t_mix)));

if reactive_transport == true
    pfd_position_z = [pfd_precip_position, pfd_age, pfd_t_mix, pfd_Cw_event, pfd_Cw_event; pfd_position_z];
else
    pfd_position_z = [pfd_precip_position, pfd_age, pfd_t_mix, pfd_Cw_event; pfd_position_z];
end


end

end
end

