%% LAST with reactive transport routine
profile on
tic
% close all;

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% selection of additional routines
plotting = false; % switches on/off the figure showing the realtime process
reactive_transport = true; % switches on/off reactive transport
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% global parameters
rho = 1000; % water density
t_end = 86400 * 7; % end time of simulation (default: 1 day = 86400 s)
dtc = 1; % initial simulation time step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% pfd parameters (in this function "prefflow_domain" the pfd parameters can be predefined and adjusted)
[mak_big, mak_mid, mak_sml, M_eventgrid, n_mak, particle_contact_grid, particle_contact_total,...
pfd_D, pfd_Diff, pfd_dim, pfd_dz, pfd_D_particle, pfd_k, pfd_m, pfd_maxV_grid, pfd_n, pfd_r,...
pfd_v, pfd_z, rate_big, rate_mid, rate_sml] = ... % output arguments
    prefflow_domain... % function
        (); % input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% soil matrix parameters
z = (0:-0.1:-1.5)'; % size of total grid
dz = abs(diff(z)); % size of one grid element
dim = length(z);
mob_fak = 0.1; % mobile particle fraction within soil matrix 

if n_mak == 0 % number of particles in the matrix dependent on the presence of macropores
   N = 1000000;
else  
   N = 2000000;
end

nclass = 800; % number of bins to subdivide the diffusivity among the different particles

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% reactive transport parameters of matrix
if reactive_transport == true
substance = 'Isoproturon';
RT_depth = -0.5; % maximum depth for reactive transport and distribution of DT50 and KF values

halflife = [3,12]; % half-life time in days, first value for topsoil and second value for subsoil
DT50 = halflife(2) * ones(dim,1);
if halflife(2) - halflife(1) ~= 0
    DT50(1:find(z==RT_depth)) = [halflife(1):(halflife(2) - halflife(1)) / (find(z==RT_depth)-1):halflife(2)]'; % depth distribution of velocity coefficient of degradation process
end
koeff = log(2) ./ DT50; %  mean velocity coefficient of degradation process

sorption_koeff = [27,3]; % Freundlich Isotherm coefficient in mL/g (taken from PPDB), first value for topsoil and second value for subsoil
KF = sorption_koeff(2) * ones(dim,1);
if sorption_koeff(2) - sorption_koeff(1) ~= 0   
    KF(1:find(z==RT_depth)) = [sorption_koeff(1):(sorption_koeff(2) - sorption_koeff(1)) / (find(z==RT_depth)-1):sorption_koeff(2)]'; % depth distribution of Freundlich coefficient of sorption process
end

beta = 0.8; % Freundlich Isotherm coefficient
rho_bulk = 1500; % soil bulk density in kg/m³
solubility = 70.2 / 1000; % solubility of substance (kg/m³), input in mg/L or g/m³


% reactive transport parameters of pfd (no depth distribution, equal
% conditions in the entire macropores)
pfd_halflife = 10; % half-life time in days in the pfd
pfd_DT50 = pfd_halflife * ones(pfd_dim,1);
pfd_koeff = log(2) ./ pfd_DT50; % mean velocity coefficient of degradation process in pfd

pfd_sorption_koeff = 10;
pfd_KF = pfd_sorption_koeff * ones(pfd_dim,1);

end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialisation script
% content: load input data; preallocating of
% arrays; conversion theta <-> psi; creation of lookup table;
% initialisation of particle measurements, positions and concentrations

init_LAST;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Start of main model process with the routines for infiltration, displacement and mixing
m_slt_surface = 0.0063; % only for simulation of IPU transport at site P4; applied IPU mass

while time < t_end 
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% initialisation of parameters v and D
v = 1 .* k; % advective velocity in each grid element of the soil matrix
D = 1 .* k ./ c; % diffusivity in each grid element of the soil matrix

v(theta < 1.1*thr(istyp)) = 0;
D(theta < 1.1*thr(istyp)) = 0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Time step and tmix
% adaptive, dynamic time stepping
%  dtc = 120;
dtc = timestep(time, t_boundary, m_surface, m, n_mak, qb_u, pfd_z, pfd_dim, pfd_v, v, D, dz, prec_time, dtc);

time = time + dtc;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tmix distribution, in each time step a new tmix value for the new
% entering particles based on this time step, ks and Dmix is calculated
[t_mix_ratio, t_mix_mean1, t_mix_SD1, t_mix_mean2, t_mix_SD2] = tmix_distribution(D_table, D, dim, ks, dtc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Infiltration routine

[bla, bla1, mass_totalsoluteinput, mass_totalwaterinput, m_surface, m_slt_surface,...
pfd_position_z, position_z_event, solute_test, solute_test2, theta] =... % output arguments
    infilt_routine... % function
        (istyp,reactive_transport, pfd_dim, t_mix_mean2,t_mix_SD2,t_mix_mean1,t_mix_SD1,t_mix_ratio, bla, bla1, Cw_eventparticles, dtc, dz, k, ks, m, mass_totalsoluteinput, mass_totalwaterinput,... % input arguments
        m_slt_surface, m_surface, n_mak, pfd_dz, pfd_k, pfd_m, pfd_position_z, pfd_r, pfd_theta, pfd_z,...
        position_z_event, psi, qb_u, rho, solute_test, solute_test2, theta, ths, time, z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Displacement routine
% based on the parameters v and D the particle displacement is firstly calculated for just the half time step (0.5*dtc) 
% and for this new time, also a new soil mositure profile is calculated from 
% which new v and D values for the second half of the time step in the corrector step can be derived

% loop pc runs two times, once for predictor and once for corrector
for pc = 1:2
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% displacement of pre-event particles in soil marix
[Cw, mtx_avg_age, position_z, theta] = ... % output arguments
    particle_track_new_mob... % function
        (D, dim, (0.5*dtc), D_table, dz, istyp, K_table, m, mob_fak, position_z, ths, v, z); % input arguments
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% displacement of event particles in soil matrix
if ~isempty(position_z_event)    
   [position_z_event, theta_event] = ... % output arguments
        particle_track_new_event... % function
            (ks, dim, (0.5*dtc), D_table, istyp, K_table, m, dz, position_z_event, prob, ths, z); % input arguments       
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% displacement of particles in pfd and drainage; only one times in the
% predictor step because it is just dependent on constant advective veloscity
if pc == 1 && ~isempty(pfd_position_z)
[pfd_particles, pfd_Cw, pfd_theta, pfd_position_z] = ... % output arguments
    particle_track_new_pfd... %function
        (n_mak, pfd_dim,pfd_dz, pfd_m, pfd_n, pfd_position_z, pfd_r, pfd_z); % input arguments
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% update of parameters psi, c and k for the corrector step
[psi, c] = psi_theta(alph, dim, istyp, n_vg, stor, theta, thr, ths); 
k = k_psi(alph, dim, istyp, ks, l_vg, n_vg, psi);

end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% reactive transport: retardation and degradation in matrix and pfd

% matrix
if reactive_transport == true
   [Cw_ret, Ct, Ct_deg, Ct0, Cm, position_z] = func_reactive_transport(dim, position_z, z, solubility, m, rho, Ct,Ct_deg, dz, theta, KF, beta, koeff, dtc); 

% pfd
    [pfd_Cw_ret, pfd_Ct, pfd_Ct_deg, pfd_Ct0, pfd_Cm, pfd_position_z] = func_pfd_reactive_transport(pfd_dim, pfd_position_z, pfd_z, solubility, pfd_m, rho, pfd_Ct,...
                                                                        pfd_Ct_deg, pfd_dz, pfd_theta, pfd_KF, beta, pfd_koeff, dtc); 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plots soil moisture and concentration conditions in soil matrix and pfd at current time step

if plotting == true
   plot_currstate(pfd_n, n_mak, pfd_m, pfd_Ct_deg, pfd_dz, pfd_Ct, pfd_Cw_ret, dz, Ct0, Cm, Ct_deg, Ct, Cw_ret, Cw, dim, istyp, pfd_Cw, pfd_dim, pfd_theta, pfd_z, theta, theta_event, ths, time, z,z_plot);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mixing between pfd and soil matrix based on the macropore depth distribution
        
[position_z_event,pfd_position_z, position_z] = ... % output arguments
    mixing_pfd_mtx... % function
        (t_mix_ratio, t_mix_SD1, t_mix_SD2, t_mix_mean1, t_mix_mean2, position_z_event,reactive_transport, dtc, k, ks, m, mak_sml, mak_mid, n_mak, particle_contact_grid,... % input arguments
        pfd_dim, pfd_dz, pfd_m, pfd_n, pfd_particles,...
        pfd_position_z, pfd_r, pfd_theta, pfd_z, position_z, psi, rate_big, rate_mid, rate_sml, rho, theta, ths, z);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Mixing of event particle with pre-event particles within the first grid elements of soil matrix
if ~isempty(position_z_event)
[position_z, position_z_event] = mixing_mtx(dim, z, reactive_transport, position_z, position_z_event, time);
end

% [val_pos, idx] = sort(position_z(:,1),'descend');
% if reactive_transport == true
%     position_z = [val_pos, position_z(idx,2), position_z(idx,3), position_z(idx,4)];
% else
%     position_z = [val_pos, position_z(idx,2), position_z(idx,3)];
% end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Updates input conditions at top of soil

ipos = find(t_boundary <= time, 1, 'last' ); % finds actual position in boundary conditions time series
qb_u = -prec_int(ipos);
Cw_eventparticles = conc(ipos);

i_time = i_time + 1;
MM(i_time) = getframe;
     
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Plots the final states of soil moisture, concentration and mass profile in soil matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% soil moisture profile
% plot_theta(D, dim, dz, theta, theta_init, theta_event, time, v, z)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% concentration profile
% plot_Cw(dz, Ct_deg, Ct, Cw_ret, concfin, Cw, Cw_init, dim, theta, theta_init, z, z_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% concentration profiles of metabolite and solute in adsorbed phase
% plot_Cm_Ct(Ct0, Cm, theta, z_plot, dim, dz)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% age profile
% plot_age(dim, mtx_avg_age, z_plot)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% mass balance and profile
plot_mass(pfd_Ct_deg, pfd_Ct, Ct, Ct_deg, Cw_ret, prec_int, conc, concfin, theta, dim, Cw, dz, theta_init, Cw_init, z, z_plot )

toc
profile viewer