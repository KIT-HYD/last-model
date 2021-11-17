%% Function for calculation of event particle displacement in the first grid elements

function [position_znew_event, theta_event] = ... % output arguments
            particle_track_new_event... % function
            (ks,dim, dt, D_table, istyp, K_table, m, dz, position_z_event, prob, ths, z) % input arguments

%% Initialisation of parameters and arrays

position_znew_event = position_z_event; % position of every event particle in soil matrix
particles = zeros(dim,1); % amount of particles in each grid element of soil matrix
theta_event = zeros(dim,1); % event soil moisture profile of soil matrix in the first grid elements
M_layer=zeros(dim,1); % total water mass in each layer

ip_0 = position_z_event(:,1) > z(1);  % correction for unphysical particle losses at top of the soil matrix
position_z_event(ip_0,1) = z(1) - position_z_event(ip_0,1); % ensures that every particle has a negative position value within the soil matrix

%% Particle displacement

   for j = 1:dim-1 
       ip = find(position_z_event(:,1) <= z(j) & position_z_event(:,1) > z(j+1)); % finds all particles within the actual grid element (between z(j) and z(j+1))
       
       if ~isempty(ip)
          
          % calculation of v and D
          v_part = ks(istyp(j)); % ks as the velocity in matrix event domain
          D_part = quantile(D_table(istyp(j)),prob); % takes the highest D value depending on the "prob-quantile" for calculating diffussivity
          
          % calculation of random walk term
          correction = 0; % correction term for random walk 
          randomwalk = 1 .* 1 .* (2*rand(length(ip),1)-ones(length(ip),1)) * (6*D_part*dt) .^ 0.5; % random walk term 
          position_znew_event(ip,1) = position_z_event(ip,1) - (((v_part+correction) * ones(length(ip),1)) * dt + randomwalk); %calculates the new positions of event water particles for the first grid elements
       
       end
       
   end

% postprocessing for realistic particle positions
ip_0 = position_znew_event(:,1) >= z(1); % correction for unphysical particle losses at top of the soil matrix due to the displacement process
position_znew_event(ip_0,1) = z(1) - position_znew_event(ip_0,1);

ip_low = position_znew_event(:,1) <= z(dim); % correction for unphysical particle losses at bottom of the soil matrix due to the displacement process
position_znew_event(ip_low,1)= z(dim-1); 

%% Calculation of water and solute masses in each grid element and generating of new soil moisture and concentration profile

for i = 1:dim-1
    ipart = find(position_znew_event(:,1) <= z(i) & position_znew_event(:,1) > z(i+1));
    particles(i) = length(ipart); % amount of particles in grid element
    
    % check oversaturation
    if particles(i) > round(ths(istyp(i))*dz(i)*1*1000/m) % if the amount of particles is higher than the saturated threshold, the amount is set to this threshold
       particles(i) = round(ths(istyp(i))*dz(i)*1*1000/m);
       
       for ii = round(ths(istyp(i))*dz(i)*1*1000/m)+1:length(ipart) % the residual particles are then moved to the next grid element below
           position_znew_event(ipart(ii),1) = position_znew_event(ipart(ii),1) - dz(i);
       end
       
    end
    
   % calculation of new event moisture profile 
   M_layer(i) = particles(i) * m; % total mass in each layer (kg)
   theta_event(i) = M_layer(i) / (dz(i)*1*1000); % soil moisture in each layer (kg/m^3)
   
%    [val_pos, idx] = sort(position_znew_event(:,1),'descend');
% position_znew_event = [val_pos, position_znew_event(idx,2), position_znew_event(idx,3), position_znew_event(idx,4), position_znew_event(idx,5)];
%    
end
 
end
