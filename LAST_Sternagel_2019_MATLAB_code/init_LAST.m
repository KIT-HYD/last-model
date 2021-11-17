%% Import dataset paths and allocate them to variables

% contains initial soil moisture state
inistate_file = 'moist_init.csv';
i_init = 1; % defines type of initial data (1= theta; 2=psi)
%contains initial concentration profile of soil matrix
iniconc_file = 'conc_init.csv';
% contains soil hydraulic parameters and depth assigment
soilpara_file = 'soilpara.dat'; 
% contains precipitation data
precip_file = 'prec_int.csv';
%contains tracer/solute concentration of precipitation water
conc_file = 'prec_conc.csv';
% contains potential evaporation data
% etp_file = 'zero_etp.dat';
% file containing observed concentration profile of soil matrix
end_conc_file = 'conc_final.csv';


%% Preallocates arrays regarding the particles of both soil matrix and pfd and others

position_z_event = []; % position of event particles entering the soil matrix
age_event = []; % age of event particles particles entering the soil matrix
Cw_event = []; % concentration of event particles entering the soil matrix

pfd_position_z = nan(1,2); % position of particles within the pfd
pfd_age = []; % age of particles within the pfd
pfd_Cw_event = []; % concentration of particles entering pfd

m_surface = 0; % water surface storage
m_slt_surface = 0; % solute surface storage
m_input = 0;

    % just for testing mass consistency
    mass_totalwaterinput = 0; % counter for amount of incoming precipitation water
    mass_totalsoluteinput = 0; % counter for amount of incoming solute mass
    
    bla = 0; % counter for water mass entering the soil matrix
    bla1 = 0; % counter for water entering the pfd
    
    solute_test = 0; % counter for solute mass entering the soil matrix
    solute_test2 = 0; % counter for solute mass entering the pfd
    
theta_event = 0;
 
%% Reads in the initial model data
    
[prec_int, prec_time] = read_precip_new(precip_file); % precipitation data (intensity and time stepping)
conc = read_conc_new(conc_file); % tracer concentration data of precipitation water
[alph, istyp, ks, l_vg, n_vg, stor, thr, ths] = soilpara(soilpara_file,dim,z); % soil data
[theta_init, psi_init] = readini_state_interp(inistate_file, i_init, z,dim); % initial theta or psi state data
Cw_init = readiniconc_state_interp(iniconc_file,z,dim); % initial tracer concentration data of soil matrix
concfin = dlmread(end_conc_file,' '); % final observed tracer concentration profile of soil matrix

%% Converts initial theta or psi state 

if i_init == 1 % calculates the initial states of psi and c from initial states of theta
    [psi_init, c_init] = psi_theta(alph, dim, istyp, n_vg, stor, theta_init, thr, ths);
else
    i_init = 2; % calculates the initial states of theta and c from initial states of psi
    [theta_init, c_init] = theta_psi(alph,dim, istyp, n_vg, psi_init, stor, thr, ths);
end


%% Creates lookup table for D and k and calculates initial hydraulic conductivity of soil matrix
    
prob = 1;
D_table = zeros(length(thr),nclass); % lookup table for D
K_table = zeros(length(thr),nclass); % lookup table for k_init
theta_table = (thr:(ths-thr)/nclass:ths); % array with soil moisture bins

%calculates for every soil moisture bin the matrix potential(psi_h), water
%capacity (c_h), diffusive coefficient (D_table) and hydraulic
%conductivity (K_table)
for jj = 1:length(thr) % create lookup table for D and k_init

    for ic=1:nclass
        [psi_h, c_h] = psi_theta(alph, dim, istyp, n_vg, stor, theta_table(ic)*ones(dim,1), thr, ths);
        k_help = k_psi(alph, dim, istyp, ks, l_vg, n_vg, psi_h);

            if c_h(1)> 0
                D_table(jj,ic) = k_help(1) / c_h(1);
            else
                D_table(jj,ic) = 0;
            end

        K_table(jj,ic) = k_help(1);
    end

end

% calculates hydraulic conductivity (initial) for every grid element of soil domain from matrix potential (psi_init)
k_init = k_psi(alph, dim, istyp, ks, l_vg, n_vg, psi_init); 

%% Initialisation of particle tracking

% sets the working parameters to their initial values
k = k_init;
c = c_init;
psi = psi_init;
theta = theta_init;
Cw = Cw_init;
pfd_Cw = pfd_Cw_initial;

% initialises particle mass (m) distribution (n), positions (position_z) and concentrations (c_particle)
[n, m] = init_particles(theta_init, dz, dim, N); 

[position_z, c_particle] = init_particles_pos(z, dim, n, Cw_init);

% sets times
time = prec_time(1); %start time (usually 0)
t_end = t_end + time;
t_boundary = prec_time; % time stepping of boundary data (time series of the precip_file)
i_time = 1;

% set initial input conditions at top soil
Cw_eventparticles = conc(1); % initial concentration of the new incoming particles
qb_u = -0.5*(prec_int(1)+prec_int(2)); % flux density of precipitation water input

%% Plot settings

z_plot = zeros((length(z)-1),1); % Changes the increments of the simulated soil grid to those of the real observed one (for a better plot fitting)
for i=1:(length(z)-1) % 
    z_plot(i) = (z(i)+z(i+1))/2;
end
z_plot=[z(1);z_plot];