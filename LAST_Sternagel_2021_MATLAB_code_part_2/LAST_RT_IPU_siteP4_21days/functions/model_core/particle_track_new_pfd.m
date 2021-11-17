%% Funtion for calculation of particle displacement in the pfd

function [pfd_particles, pfd_Cw, pfd_theta, pfd_position_znew] = ... % output arguments
            particle_track_new_pfd... %function
            (n_mak, pfd_dim,... % input arguments
            pfd_dz, pfd_m, pfd_n, pfd_position_z, pfd_r, pfd_z)

%% Initialisation of parameters and arrays

pfd_position_znew = pfd_position_z; % position of every particle in pfd
pfd_particles = zeros(pfd_dim,1); % amount of pfd_particles in each grid element of pfd
pfd_theta = zeros(pfd_dim,1); % soil moisture profile of pfd
pfd_Cw = zeros(pfd_dim,1); % concentration profile of pfd
pfd_age = zeros(pfd_dim,1); % average age progile of pfd (only for testin, yet)
M_layer = zeros(pfd_dim,1); % total mass in each layer

ip_0 = find(pfd_position_z(:,1) > pfd_z(1));  % correction for unphysical particle losses at top of the pfd
pfd_position_z(ip_0,1) = pfd_z(1) - pfd_position_z(ip_0,1); % ensures that every particle has a negative position value within the pfd

%% Particle displacement

   for j = 1:pfd_dim-1 
       ip = find(pfd_position_z(:,1) <= pfd_z(j) & pfd_position_z(:,1) > pfd_z(j+1)); % finds all pfd_particles within the actual grid element (between z(j) and z(j+1))
  
       if ~isempty(ip)   
            pfd_position_znew(ip,1) = pfd_z(pfd_dim);     
       end
       
   end

%% Calculation of water and solute masses in each grid element and generating of new soil moisture and concentration profile

for i = pfd_dim-1:-1:1
    ipart = find(pfd_position_znew(:,1) < pfd_z(i) & pfd_position_znew(:,1) >= pfd_z(i+1));
    pfd_particles(i)=length(ipart); % amount of pfd_particles in grid element
    
    % check oversaturation
    if pfd_particles(i) > pfd_n*n_mak % if the number of pfd_particles is higher than the saturated threshold the number is set to this threshold
       pfd_particles(i) = round(pfd_n(i)*n_mak);
       
       for ii= ipart(1):(ipart(end) - round(pfd_n(i)*n_mak)) % redistribution of residual pfd_particles out of the last grid element to the residual unsaturated grid elements above
           pfd_position_znew(ii,1)= pfd_position_znew(ii,1) + pfd_dz(i); 
       end
       
    end
    
    ip_0 = find(pfd_position_znew(:,1) > pfd_z(1));  % correction for unphysical particle losses at top of the pfd
    pfd_position_znew(ip_0,1) = pfd_z(1);

    blabla = find(pfd_position_znew(:,1) == pfd_z(i+1));
    pfd_position_znew(blabla,1) = pfd_z(i+1) + (pfd_z(i)-pfd_z(i+1)) .* rand(length(blabla),1); % every particle gets randomly an own position in the respective grid element (just for style)
    
    % calculation of new profiles
    M_layer(i) = pfd_particles(i) * pfd_m; % total mass in each grid element
    pfd_theta(i) = round(M_layer(i)/((pfd_dz(i)*(pi*pfd_r^2)*1000)*n_mak),3); % new soil moisture in each grid element

    ipart2 = find(pfd_position_znew(:,1) < pfd_z(i) & pfd_position_znew(:,1) >= pfd_z(i+1));
    pfd_age(i) = sum(pfd_position_znew(ipart2,2))/pfd_particles(i); % new average age in each grid element
    pfd_Cw(i) = sum(pfd_position_znew(ipart2,4))/pfd_particles(i); % new concentration in each grid element
    
        if isnan(pfd_Cw(i))
           pfd_Cw(i) = 0; 
        end 

        if isnan(pfd_age(i))
           pfd_age(i) = 0;
        end 
    
    pfd_position_znew(ipart2,4) = pfd_Cw(i) * ones(length(ipart2),1); % new concentration of every single particle in each grid element
end

end
