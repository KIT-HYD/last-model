%% Function for defining and calculating the parameters of the pfd
% in the first two initialisation blocks the most important parameters can
% be predefinied and adjusted

function [mak_big, mak_mid, mak_sml, M_eventgrid, n_mak, particle_contact_grid, particle_contact_total,...
          pfd_D, pfd_Diff, pfd_dim, pfd_dz, pfd_D_particle, pfd_k, pfd_m, pfd_maxV_grid, pfd_n, pfd_r,...
          pfd_v, pfd_z, rate_big, rate_mid, rate_sml] = ... % output arguments
            prefflow_domain... % function
            () % input arguments

%assumption: preferential flow domain with certain number of macropores
%shaped like a straight circular cylinder.
%Water particles are spherical shaped in a cubic storage

%% Initialisation of pfd and macropore geometry

pfd_z = (0:-0.05:-0.20)'; % maximum length of a macropore
pfd_D = (3/(10*100)); % diameter of macropore (mm)
n_mak = 68; % amount of macropores within the domain

pfd_r = pfd_D/2; % radius of a macropore
pfd_dz = abs(diff(pfd_z)); % length of a macropore grid element
pfd_dim = length(pfd_z);
pfd_N = 10000; % amount of particles within one macropore
rho = 1000; % density water 

%% Initialisation of macropore depth distribution

mak_big = -0.20; % depth distribution of macropores within pfd
mak_mid = -0.0;
mak_sml = -0.0;

rate_big = 1.0; % distribution factors for the proportion of macropores and diffusive mixing masses
rate_mid = 0.0; % when adjusting the depth of the biggest macropores, please do it by adjusting "pfd_z" above
rate_sml = 0.0;

%% Calculation of further pfd/macropore parameters

pfd_maxV_grid = pi * (pfd_r^2) * abs(pfd_dz); % volume of grid element
V_gesamt = sum(pfd_maxV_grid); % total volume of a macropore

M_eventgrid = rho * pfd_maxV_grid; % water mass fitting in a grid element
m_gesamt = sum(M_eventgrid); % total water mass fitting in a macropore 
pfd_m = m_gesamt / pfd_N; % mass of one particle

V_particle = pfd_m/rho; % volume of a particle 
pfd_D_particle = nthroot((V_particle/(pi/6)),3); % diameter of a particle
r_particle = pfd_D_particle / 2; % radius of a macropore

pfd_n = M_eventgrid / pfd_m; % amount of particles fitting into a grid element
n_gesamt = sum(pfd_n); % amount of particles fitting into a macropore

U_eventgrid = 2 * pi * pfd_r; % circumference of a grid element

particle_contact = floor(U_eventgrid/pfd_D_particle); % amount of particles fitting with their diameter next to each other in a row 
anzahl_reihen = floor(abs(pfd_dz(1))/pfd_D_particle); % amount of rows fitting over each other 
particle_contact_grid = particle_contact * anzahl_reihen; % amount of particles having contact to lateral surface of grid element
particle_contact_total = particle_contact_grid * length(pfd_dz); % total amount of particles having contact to lateral surface of a macropore

if particle_contact_total > pfd_N % ensures that the amount of contact particles is not higher than the total amount of particles within the pfd
   particle_contact_total = pfd_N;
   particle_contact_grid = particle_contact_total / length(pfd_dz);
end

pfd_qmak = 2884.2 *(pfd_r^2); % flux density of a macropore
pfd_v = pfd_qmak; % advective velocity of a particle 
pfd_k = pfd_v; % hydraulic conductivity 
pfd_Diff = 0; % no diffusion within macropore

end

