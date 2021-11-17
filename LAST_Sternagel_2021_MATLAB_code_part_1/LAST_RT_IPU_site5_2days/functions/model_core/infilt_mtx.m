%% Function to describe the infiltraion into the soil matrix

function [bla, m_slt_surface, m_surface, position_z_event, solute_test, theta, theta_count] = ... % output arguments
            infilt_mtx... % function
            (t_mix_mean2,t_mix_SD2,t_mix_mean1,t_mix_SD1,t_mix_ratio, bla, Cw_eventparticles, dz, m, m_input1, m_slt_surface, m_surface, position_z_event,...
            rho, solute_test, theta, time, z) % input arguments


  
precip_particles = floor((m_input1)/m); % amount of event particles to be injected, converts the mass input to number of particles
bla = bla + (precip_particles*m); % counter for amount of water totally infiltrating the soil matrix
m_surface = (m_surface - m_input1) + (m_input1 - (precip_particles*m)) ; % updates water surface storage
m_slt_surface = (m_slt_surface - (Cw_eventparticles*(m_input1/rho))) + (Cw_eventparticles*((m_input1 - (precip_particles*m))/rho)); % updates solute mass surface storage
solute_test = solute_test + (Cw_eventparticles*((precip_particles*m)/rho)); % counter for amount of solutes totally infiltrating the soil matrix within the first scenario
precip_position = z(1) * ones(precip_particles,1); % in every time step the new infiltrating particles get the position z(1)
theta_count = 1; % sets scenario counter to 1, thus the routine will not run through the subsequent scenarios
theta(1) = theta(1) + precip_particles*m/(rho*dz(1)*1); % updates the soil moisture of the first grid element
age_event = time*ones(precip_particles,1); % matrix with the age of the infiltrating event particles, every new particle gets the value of its timestep when entering the domain
 

t_mix_event = [randn(floor(precip_particles*t_mix_ratio),1) * t_mix_SD1 + t_mix_mean1;randn(ceil(precip_particles*(1-t_mix_ratio)),1) * t_mix_SD2 + t_mix_mean2];
t_mix_event = t_mix_event(randperm(numel(t_mix_event)));


     %adds the number of event particles entering the soil matrix with
     %their concentration value to the Cw_event matrix
     if Cw_eventparticles == 0
        Cw_event = zeros(precip_particles,1);      
     else 
        Cw_event = Cw_eventparticles*ones(precip_particles,1);
     end
     
position_z_event = [precip_position, age_event, t_mix_event, Cw_event, Cw_event; position_z_event];





end

