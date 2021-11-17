%% Description:
% This code comprises routines considering solute transport under
% preferential flow conditions and mixing between macropores and soil
% matrix and therewith extending the original particle based
% Lagrangian Model of Zehe and Jackisch (2016). This extension was done
% within the scope of the master thesis of Alexander Sternagel (2018).
% Further information on the theoretical background and the development can
% be obtained from these two references. The functionality and handling of
% this code/model are explained in the different readme-files within the
% folder structure.
% 
% Zehe, E. and Jackisch, C.: A Lagrangian model for soil water dynamics during rainfall-driven conditions,
% Hydrol. Earth Syst. Sci., 20, 3511–3526, doi:10.5194/hess-20-3511-2016, 2016
% 
% Sternagel, A.: Modelling of Solute Transport and Travel Time
% Distributions Using a Particle-Based Lagrangian Model, Master Thesis, 2018

profile on
tic
close all;

%% Paths to the different folders containing all the data, functions and scripts

addpath('source\boundarycon')
addpath('source\initial_states')
addpath('source\soils')

addpath('functions\initial')
addpath('functions\load_model_data')
addpath('functions\model_core')
addpath('functions\plot_functions')
addpath('functions\soil_functions')


%% Initialisation: Predefines global, soil matrix and pfd parameters and loads the initialisation script
% here the gloabal and soil matrix parameters can be predefined and
% adjusted

% global parameters
rho = 1000; % water density
t_end = 86400; % end time of simulation (default: 1 day = 86400 s)
dtc = 120; % simulation time step

% dtc = ceil(abs(pfd_z(pfd_dim))/pfd_v); % alternative: dynamic time step dependent on advective velocity in pfd     
% if dtc < 20 % ensures a not to small time stepping 
%    dtc = 20; 
% end

% soil matrix parameters
z = (0:-0.1:-1.5)'; % size of total grid
dz = abs(diff(z)); % size of one grid element
dim = length(z);
mob_fak = 0.1; % mobile particle fraction within soil matrix 
N = 1000000; % total number of particles within soil matrix
nclass = 800; % number of bins to subdivide the diffusivity among the different particles
t_mix = 1; % predefining the mixing time between event and pre-event particles within soil matrix

% pfd parameters (in this function "prefflow_domain" the pfd parameters can be predefined and adjusted)
[mak_big, mak_mid, mak_sml, M_eventgrid, n_mak, particle_contact_grid, particle_contact_total, pfd_Cw_initial,...
pfd_D, pfd_Diff, pfd_dim, pfd_dz, pfd_D_particle, pfd_k, pfd_m, pfd_maxV_grid, pfd_n, pfd_r, pfd_theta,...
pfd_v, pfd_z, rate_big, rate_mid, rate_sml] = ... % output arguments
    prefflow_domain... % function
        (); % input arguments
    
    % Initialisation script
    % content: load input data; preallocating of
    % arrays; conversion theta <-> psi; creation of lookup table;
    % initialisation of particle measurements, positions and concentrations

    init_LAST;
    
drainage = false; % switches on/off the drainage at the bottom of the macropores
plotting = false; % switches on/off the figure showing the realtime process

%% Start of main model process with the routines for infiltration, displacement and mixing

while time < t_end 

time = time + dtc;
 
%% Infiltration routine

[age_event, bla, bla1, Cw_event, mass_totalsoluteinput, mass_totalwaterinput, m_input, m_surface, m_slt_surface, pfd_age, pfd_Cw_event,...
pfd_position_z, position_z_event, solute_test, solute_test2, theta] =... % output arguments
    infilt_routine... % function
        (age_event, bla, bla1, Cw_event, Cw_eventparticles, dtc, dz, k, ks, m, mass_totalsoluteinput, mass_totalwaterinput, m_input,... % input arguments
        m_slt_surface, m_surface, n_mak, pfd_age, pfd_Cw_event, pfd_dz, pfd_k, pfd_m, pfd_position_z, pfd_r, pfd_theta, pfd_z,...
        position_z_event, psi_init, qb_u, rho, solute_test, solute_test2, theta, ths, time, z);

%% Displacement routine
% based on the parameters v and D the particle displacement is firstly calculated for just the half time step (0.5*dtc) 
% and for this new time, also a new soil mositure profile is calculated from 
% which new v and D values for the second half of the time step in the corrector step can be derived

% loop pc runs two times, once for predictor and once for corrector
for pc = 1:2
% initialisation of parameters v and D
v = 1.*k; % advective velocity in each grid element of the soil matrix
D = 1.*k./(c); % diffusivity in each grid element of the soil matrix

    ipres = find(theta < 1.1*thr(istyp)); % finds all grid elements with a soil moisture near to thr (almost dry soil)
    if length(ipres) >= 1 % sets the velocity and diffusivity in this grid elements to 0-> no flux!
        v(ipres) = 0;
        D(ipres) = 0;
    end

% displacement of pre-event particles in soil marix
[c_particle, Cw, position_z, theta] = ... % output arguments
    particle_track_new_mob... % function
        (c_particle, D, dim, (0.5*dtc), D_table, dz, istyp, K_table, m, mob_fak, position_z, ths, v, z); % input arguments
        
% displacement of event particles in soil matrix    
if ~isempty(position_z_event) 
    
   [position_z_event, theta_event] = ... % output arguments
        particle_track_new_event... % function
            (dim, (0.5*dtc), D_table, istyp, K_table, m, dz, position_z_event, prob, ths, z); % input arguments
        
end

% displacement of particles in pfd and drainage; only one times in the
% predictor step because it is just dependent on constant advective veloscity
if pc == 1
[c_particle, mtx_particles_mix0, pfd_particles, pfd_particle_output, pfd_age, pfd_Cw, pfd_Cw_event, pfd_theta, pfd_position_z, position_z] = ... % output arguments
    particle_track_new_pfd... %function
        (c_particle, drainage, dtc, k, m, mak_big, mak_mid, mak_sml, n_mak, pfd_age, pfd_Cw_event, pfd_Diff, pfd_dim,... % input arguments
        pfd_dz, pfd_m, pfd_n, pfd_position_z, pfd_r, pfd_v, pfd_z, position_z, rate_big, rate_mid, rate_sml, z);
end

% update of parameters psi, c and k for the corrector step
[psi, c] = psi_theta(alph, dim, istyp, n_vg, stor, theta, thr, ths); 
k = k_psi(alph, dim, istyp, ks, l_vg, n_vg, psi);

end
              
%% Plots soil moisture and concentration conditions in soil matrix and pfd at current time step

if plotting == true
   plot_currstate(Cw, dim, istyp, pfd_Cw, pfd_dim, pfd_theta, pfd_z, theta, theta_event, ths, time, z),
end
        
%% Mixing of event particle with pre-event particles within the first grid elements of soil matrix

[c_particle, Cw_event, position_z, position_z_event] = mixing_mtx(age_event, c_particle, Cw_event, position_z, position_z_event, time, t_mix);

%% Mixing between pfd and soil matrix based on the macropore depth distribution
        
[c_particle, pfd_age, pfd_Cw_event, pfd_position_z, position_z] = ... % output arguments
    mixing_pfd_mtx... % function
        (c_particle, dtc, k, ks, m, mak_sml, mak_mid, mtx_particles_mix0, n_mak, particle_contact_grid,... % input arguments
        pfd_age, pfd_Cw_event, pfd_dim, pfd_dz, pfd_m, pfd_n, pfd_particles, pfd_particle_output,...
        pfd_position_z, pfd_r, pfd_theta, pfd_z, position_z, psi, rate_big, rate_mid, rate_sml, rho, theta, ths, z);       

%% Updates input conditions at top of soil

ipos = find(t_boundary <= time, 1, 'last' ); % finds actual position in boundary conditions time series
qb_u = -0.5 * (prec_int(ipos)+prec_int(ipos+1)); % updates flux density of precipitation water input
Cw_eventparticles = conc(ipos); % updates the concentration of precipitation water input
i_time = i_time+1;
MM(i_time) = getframe;
     
end 
     
%% Plots the final states of soil moisture, concentration and mass profile in soil matrix

% soil moisture profile
plot_theta(D, dim, dz, theta, theta_init, theta_event, time, v, z)

% concentration profile
plot_Cw(concfin, Cw, Cw_init, dim, theta, theta_init, z, z_plot)

% mass balance and profile
plot_mass(prec_int, conc, concfin, theta, dim, Cw, dz, theta_init, Cw_init, z, z_plot)

toc
profile viewer