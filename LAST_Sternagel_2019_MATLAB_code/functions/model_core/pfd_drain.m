%% Function for calculation of a drainage at the bottom of the macropores

function [c_particle, pfd_particles, pfd_particle_output, pfd_age, pfd_Cw, pfd_Cw_event, pfd_theta, pfd_position_znew, position_z] = ... % output arguments
            pfd_drain... % function
             (c_particle, dt, k, m, mak_big, mak_mid, mak_sml, M_layer, n_mak, pfd_age, pfd_Cw_event, pfd_dim, pfd_dz, pfd_m,... % input arguments
             pfd_position_znew, pfd_r, pfd_theta, pfd_v, pfd_z, position_z, rate_big, rate_mid, rate_sml, z)

%% Initialisation of parameters and arrays

pfd_Cw_event2 = zeros(pfd_dim,1); % assistant parameter for calculation of new event concentration profile of pfd
pfd_particles = zeros(pfd_dim,1); % amount of particles in each grid element of pfd
avg_age = zeros(pfd_dim,1); % average age progile of pfd (only for testing, yet)

%% Calculation of masses and particles leaving the pfd at the bottom

kp_indx = find(pfd_z(pfd_dim-1) > z); % finds respective depth of soil matrix corresponding to the depth of the pfd
pfd_t_lb = (pfd_r*2) * (sqrt(1/pfd_v * 1/k(kp_indx(1)))); % calculates mixing time for diffusion from pfd into soil matrix
m_output = (M_layer(pfd_dim-1)/pfd_t_lb) * dt; % calculates leaving mass

if m_output > M_layer(pfd_dim-1) % ensures only the maximum possible mass is leaving
   m_output = M_layer(pfd_dim-1);
end

pfd_particle_output = floor(m_output/pfd_m); % converts the mass into integer number of particles

%% Deletion of the particles out of all pfd arrays and updating of soil moisture

zzz = find(pfd_position_znew(:,1) < pfd_z(pfd_dim-1) & pfd_position_znew(:,1) >= pfd_z(pfd_dim));

pfd_theta(pfd_dim-1) = pfd_theta(pfd_dim-1) - round(m_output/((pfd_dz(pfd_dim-1)*(pi*pfd_r^2)*1000)*n_mak),3);
pfd_position_znew((zzz(end)-(pfd_particle_output-1):zzz(end)),:) = [];
pfd_age(((zzz(end)-(pfd_particle_output-1):zzz(end)))) = [];
pfd_conc_output = sum(pfd_Cw_event((zzz(end)-(pfd_particle_output-1):zzz(end))))/pfd_particle_output; % solute concentration of leaving particles
pfd_Cw_event((zzz(end)-(pfd_particle_output-1):zzz(end))) = [];

%% Assumption: If the next grid element above is then more saturated it comes immediately to a flux of particles
%% from this saturated grid element to the now unsaturated last grid element

for i = pfd_dim-1:-1:2 

    if pfd_theta(i-1) > pfd_theta(i)
       pfd_theta_diff = pfd_theta(i-1) - pfd_theta(i);
       pfd_theta(i) = pfd_theta(i) + pfd_theta_diff; 
       pfd_theta(i-1) = pfd_theta(i-1) - pfd_theta_diff;

       ipart = find(pfd_position_znew(:,1) < pfd_z(i-1) & pfd_position_znew(:,1) >= pfd_z(i));
       pfd_position_znew(ipart(end)-(pfd_particle_output-1):ipart(end)) = pfd_position_znew(ipart(end)-(pfd_particle_output-1):ipart(end)) - pfd_dz(i);
    end

end    

%% Updating of concentration and age profiles of pfd after leaving and redistribution of particles

for iii = 1:pfd_dim-1 
    ipart3 = find(pfd_position_znew(:,1) < pfd_z(iii) & pfd_position_znew(:,1) >= pfd_z(iii+1));
    pfd_particles(iii) = length(ipart3);
    pfd_Cw_event2(iii) = sum(pfd_Cw_event(ipart3));
    pfd_Cw_event(ipart3) = (pfd_Cw_event2(iii)./pfd_particles(iii))*ones(length(ipart3),1);
    avg_age(iii) = sum(pfd_position_znew(ipart3,2))/pfd_particles(iii);
end

pfd_Cw = pfd_Cw_event2 ./ pfd_particles;
pfd_Cw(isnan(pfd_Cw)) = 0;
avg_age(isnan(avg_age)) = 0;

%% Conversion of leaving pfd particles into matrix pfd_particles and distributive entering of the soil matrix in the certain depths of the three macropore sizes

% finds respective depths of soil matrix corresponding to the bottom of the
% big, medium and small macropores
indx0_big = find(pfd_z(pfd_z==mak_big) == round(z,2)); 
indx0_mid = find(pfd_z(pfd_z==mak_mid) == round(z,2)); 
indx0_sml = find(round(pfd_z(round(pfd_z,2)==mak_sml),2) == round(z,2));

mtx_particles_mix0 = floor((pfd_particle_output*pfd_m)/m); % conversion of pfd particles into soil matrix pfd_particles

% new positions of the entering particles in the three depths of the soil matrix
pos0_big = z(indx0_big(1)) + (z(indx0_big(1)-1)-z(indx0_big(1))) .* rand(round(mtx_particles_mix0*rate_big),1); 
pos0_mid = z(indx0_mid(1)) + (z(indx0_mid(1)-1)-z(indx0_mid(1))) .* rand(round(mtx_particles_mix0*rate_mid),1);
pos0_sml = z(indx0_sml(1)) + (z(indx0_sml(1)-1)-z(indx0_sml(1))) .* rand(round(mtx_particles_mix0*rate_sml),1);

position_z = sort([pos0_big;pos0_mid;pos0_sml;position_z],'descend'); % adding the new particle positions to the position array of the soil matrix

% adding the concentrations of the new particles to the three respective depths of the particle
% concentration array of the soil matrix
down_big = find(position_z <= z(indx0_big(1)-1)); 
up_big = c_particle(1:down_big-1);
c_particle(1:down_big-1) = [];
c_particle = [pfd_conc_output*ones(length(pos0_big),1);c_particle];
c_particle = [up_big;c_particle];

down_mid = find(position_z <= z(indx0_mid(1)-1)); 
up_mid = c_particle(1:down_mid-1);
c_particle(1:down_mid-1) = [];
c_particle = [pfd_conc_output*ones(length(pos0_mid),1);c_particle];
c_particle = [up_mid;c_particle];

down_sml = find(position_z <= z(indx0_sml(1)-1)); 
up_sml = c_particle(1:down_sml-1);
c_particle(1:down_sml-1) = [];
c_particle = [pfd_conc_output*ones(length(pos0_sml),1);c_particle];
c_particle = [up_sml;c_particle];
   
end

