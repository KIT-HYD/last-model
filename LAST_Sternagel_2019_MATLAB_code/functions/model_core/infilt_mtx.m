%% Function to describe the infiltraion into the soil matrix

function [age_event, bla, Cw_event, m_input, m_slt_surface, m_surface, position_z_event, solute_test, theta, theta_count] = ... % output arguments
            infilt_mtx... % function
            (age_event, bla, Cw_event, Cw_eventparticles, dz, m, m_input1, m_slt_surface, m_surface, position_z_event,...
            rho, solute_test, theta, time, z) % input arguments

if m_input1 > m
  
 precip_particles = floor((m_input1)/m); % amount of event particles to be injected, converts the mass input to number of particles
 bla = bla + (precip_particles*m); % counter for amount of water totally infiltrating the soil matrix
 m_surface = m_surface - m_input1; % updates water surface storage
 m_slt_surface = m_slt_surface - (Cw_eventparticles*m_input1); % updates solute mass surface storage
 solute_test = solute_test + (Cw_eventparticles*m_input1); % counter for amount of solutes totally infiltrating the soil matrix within the first scenario
 m_input = (m_input1) - (precip_particles*m); % calculates the remaining water particles, which where not considered by the "floor" function, and converts them back into a mass (in mm) to add it in the next time step. So no water input gets lost
 precip_position = z(1) * ones(precip_particles,1); % in every time step the new infiltrating particles get the position z(1)
 theta_count = 1; % sets scenario counter to 1, thus the routine will not run through the subsequent scenarios
 theta(1) = theta(1) + precip_particles*m/(rho*dz(1)*1); % updates the soil moisture of the first grid element
 [position_z_event] = [precip_position;position_z_event]; % matrix with the positions of the infiltrating event particles, initial position for every particle is always 0 and then they were displaced in every time step 
 [age_event] = [time*ones(precip_particles,1);age_event]; % matrix with the age of the infiltrating event particles, every new particle gets the value of its timestep when entering the domain

     %adds the number of event particles entering the soil matrix with
     %their concentration value to the Cw_event matrix
     if Cw_eventparticles == 0
        [Cw_event] = [zeros(length(position_z_event),1); Cw_event];
     else 
        [Cw_event] = [Cw_eventparticles*ones(length(position_z_event),1); Cw_event];
     end

 end

end

