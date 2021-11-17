%% Funtion for calculation of particle displacement in the pfd

function [c_particle, mtx_particles_mix0, pfd_particles, pfd_particle_output, pfd_age, pfd_Cw, pfd_Cw_event, pfd_theta, pfd_position_znew, position_z] = ... % output arguments
            particle_track_new_pfd... %function
            (c_particle, drainage, dt, k, m, mak_big, mak_mid, mak_sml, n_mak, pfd_age, pfd_Cw_event, pfd_Diff, pfd_dim,... % input arguments
            pfd_dz, pfd_m, pfd_n, pfd_position_z, pfd_r, pfd_v, pfd_z, position_z, rate_big, rate_mid, rate_sml, z)

%% Initialisation of parameters and arrays

pfd_position_znew = pfd_position_z; % position of every particle in pfd
pfd_particles = zeros(pfd_dim,1); % amount of pfd_particles in each grid element of pfd
pfd_theta = zeros(pfd_dim,1); % soil moisture profile of pfd
pfd_Cw = zeros(pfd_dim,1); % concentration profile of pfd
avg_age = zeros(pfd_dim,1); % average age progile of pfd (only for testin, yet)
M_layer=zeros(pfd_dim,1); % total mass in each layer

ip_0=find(pfd_position_z(:,1) > pfd_z(1));  % correction for unphysical particle losses at top of the pfd
pfd_position_z(ip_0)=pfd_z(1)-pfd_position_z(ip_0); % ensures that every particle has a negative position value within the pfd

%% Particle displacement

   for j = 1:pfd_dim-1 
       ip = find(pfd_position_z(:,1) <= pfd_z(j) & pfd_position_z(:,1) > pfd_z(j+1)); % finds all pfd_particles within the actual grid element (between z(j) and z(j+1))
  
       if ~isempty(ip)
          correction = 0; % correction term for random walk (no correction because no spatial variable diffussion-> all event pfd_particles travel with the same v and D through the preferential domain)
          randomwalk = 1 .* 1 .* (2*rand(length(ip),1)-ones(length(ip),1)) * (6*pfd_Diff*dt) .^ 0.5; % random walk term
          pfd_position_znew(ip) = pfd_position_z(ip) - (pfd_v+correction) * ones(length(ip),1) * dt + randomwalk; % calculates the new positions of pfd_particles
          pfd_position_znew(ip(pfd_position_znew(ip)<pfd_z(pfd_dim))) = pfd_z(pfd_dim); % ensures every particle is in the last grid element
       end
       
   end

% postprocessing for realistic particle positions
ip_0 = find(pfd_position_znew(:,1) >= pfd_z(1)); % correction for unphysical particle losses at bottom of the pfd due to the displacement process
pfd_position_znew(ip_0) = pfd_z(1) - pfd_position_znew(ip_0);

%% Calculation of water and solute masses in each grid element and generating of new soil moisture and concentration profile

for i = pfd_dim-1:-1:1
    ipart = find(pfd_position_znew(:,1) < pfd_z(i) & pfd_position_znew(:,1) >= pfd_z(i+1));
    pfd_particles(i)=length(ipart); % amount of pfd_particles in grid element
    
    % check oversaturation
    if pfd_particles(i) > pfd_n*n_mak % if the number of pfd_particles is higher than the saturated threshold the number is set to this threshold
       pfd_particles(i) = round(pfd_n(i)*n_mak);
       
       for ii= ipart(1):(ipart(end) - round(pfd_n(i)*n_mak)) % redistribution of residual pfd_particles out of the last grid element to the residual unsaturated grid elements above
           pfd_position_znew(ii)= pfd_position_znew(ii) + pfd_dz(i); 
       end
       
    end
    
    pfd_position_znew(ipart(pfd_position_znew(ipart)>pfd_z(1))) = pfd_z(1); % again a correction for unwanted pfd_particles losses at the top of the pfd
    blabla = find(pfd_position_znew(:,1) == pfd_z(i+1));
    pfd_position_znew(blabla) = pfd_z(i+1) + (pfd_z(i)-pfd_z(i+1)) .* rand(length(blabla),1); % very particle gets randomly an own position in the respective grid element (just for style)
    
    % calculation of new profiles
    M_layer(i) = pfd_particles(i) * pfd_m; % total mass in each grid element
    pfd_theta(i) = round(M_layer(i)/((pfd_dz(i)*(pi*pfd_r^2)*1000)*n_mak),3); % new soil moisture in each grid element

    ipart2 = find(pfd_position_znew(:,1) < pfd_z(i) & pfd_position_znew(:,1) >= pfd_z(i+1));
    avg_age(i) = sum(pfd_position_znew(ipart2,2))/pfd_particles(i); % new average age in each grid element
    pfd_Cw(i)=sum(pfd_Cw_event(ipart2))/pfd_particles(i); % new concentration in each grid element
    
        if isnan(pfd_Cw(i))
           pfd_Cw(i) = 0; 
        end 

        if isnan(avg_age(i))
           avg_age(i) = 0;
        end 
    
    pfd_Cw_event(ipart2)=pfd_Cw(i)*ones(length(ipart2),1); % new concentration of every single particle in each grid element

end

%% Calculation of drainage at the bottom of pfd (optional)

    if drainage == true

    [c_particle, pfd_particles, pfd_particle_output, pfd_age, pfd_Cw, pfd_Cw_event, pfd_theta, pfd_position_znew, position_z] = ... % output arguments
        pfd_drain... % function
        (c_particle, dt, k, m, mak_big, mak_mid, mak_sml, M_layer, n_mak, pfd_age, pfd_Cw_event, pfd_dim,...
        pfd_dz, pfd_m, pfd_position_znew, pfd_r, pfd_theta, pfd_v, pfd_z, position_z, rate_big, rate_mid, rate_sml, z); % input arguments
    else
        pfd_particle_output = 0;
        mtx_particles_mix0 = 0;

    end


end
