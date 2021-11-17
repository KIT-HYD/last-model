%% DIPMI approach for simulating experiment of Bowers et al. (2020)

%% initialisation
% soil parameters
soilpara_file = 'soilpara.dat'; 
[alph, ks, l_vg, n_vg, stor, thr, ths] = soilpara(soilpara_file); 
m_vg = 1 - 1/n_vg;
bulk_density = 1200;
N_particles = 100000; % number of particles
N_bins = (1:200); % number of bins
%--------------------------------------------------------------------------
% soil water content bins
theta_bins = [ths:-(ths-thr) / length(N_bins):thr+(ths-thr) / length(N_bins)];
%--------------------------------------------------------------------------
% particle array initialisation
water_mass_total = (ths * (0.03/ bulk_density)) * 1000; % total water mass (kg) at saturation in one sample of 30 g
water_mass_particle = water_mass_total / N_particles; % water mass per particle

particle_array = zeros(N_particles,3); % pre-allocating particle array
%--------------------------------------------------------------------------
% particle distribution over bins (equal) in saturation 
bin_particles = (N_particles / length(N_bins)) * ones(1,length(N_bins));
bin_cumsum_particles = cumsum(bin_particles);

for i = length(bin_particles):-1:1
    particle_array(bin_cumsum_particles(i)-bin_particles(i)+1:bin_cumsum_particles(i)) = N_bins(i);  
end
%--------------------------------------------------------------------------
% assignment of isotopic signature/concentration to particles
% Assumption: 30 g soil with a total of approx. 1.7 mL water (20/350 * 30). 1.7 mL
% = g based on the total possible amount of water of 10.7 g at full saturation
% corresponds to approx. 16.6% relative saturation. Converted to volumetric
% The water content is theta = 0.122 after initial moistening with
% light water (signature 1). This theta corresponds to a bin occupancy
% from 200 to 168. The particles in these bins are given a signature for
% light water (-99 H and -12.3 O for upper boundary of SD; -79 H and -9.3 O for lower boundary of SD)...

for i = 168:200
    particle_array(particle_array == i,2) = -99;
    particle_array(particle_array == i,3) = -12.3;
    
%     particle_array(particle_array == i,2) = -79;
%     particle_array(particle_array == i,3) = -9.3;
end

% ...the rest to full saturation is done with heavy water (-48 H and -7.8 O for upper boundary of SD; -46 H and -7.2 O for lower boundary of SD)
for i = 1:167
    particle_array(particle_array == i,2) = -48;
    particle_array(particle_array == i,3) = -7.8;
    
%     particle_array(particle_array == i,2) = -46;
%     particle_array(particle_array == i,3) = -7.2;
end
%--------------------------------------------------------------------------
% young-laplace equation to calculate radii of different pores and to determine characteristic flow length LD
psi_bins_m = psi_theta(alph, n_vg, theta_bins, thr, ths); % matric potentials (m) in bins
sigma = 0.073; % surface tension water kg/s²
rho = 1000; % density water kg/m³
g = 9.81; % earth accelaration m/s²

%young-laplace equation
r_bins_m = -(2 * sigma) ./ (psi_bins_m * rho * g); % pore radii in m
r_bins_ym = r_bins_m * 10^6;
y1 = flip(cumsum(r_bins_ym(2:end),'reverse'));
r_bins_ym(1) = interp1([1:length(N_bins)-1],y1,[length(N_bins)],'pchip','extrap') - y1(end);
ceil(trapz(r_bins_ym));

x = [round(ceil(trapz(r_bins_ym)),-3)];
%--------------------------------------------------------------------------
% defining edges/locations of bins in LD
bin_edges = [x:-(x / length(N_bins)):0];

particle_array = [zeros(N_particles,1),particle_array];

for iii = 1:200
    particle_array(particle_array(:,2) == iii,1) = bin_edges(iii+1) + (bin_edges(iii) - bin_edges(iii+1)) * rand(sum(particle_array(:,2) == iii),1);
end
%--------------------------------------------------------------------------
% pore-size-dependent diffusivity in bins/tension areas
porosity = ths;
porosity_factor = (theta_bins - thr) ./ porosity;
D_bins = (2.299 * 10^-9 * porosity_factor);
%--------------------------------------------------------------------------
%% start diffusive mixing among bins
% diffusive mixing (model core)
t_end = 3600*24*7;
time = 0;
dt = 600;
result_array = zeros(2,length(N_bins));
C_bins = zeros(2,length(N_bins)); % isotope concentrations in bins
%--------------------------------------------------------------------------
% C_bins at t = 0
for iii = 1:200
    
    ip = (particle_array(:,2) == iii);
    
    % mean conc. H in bins
    C_bins(1,iii) = mean(particle_array(ip,3));
    
    % mean conc. O in bins
    C_bins(2,iii) = mean(particle_array(ip,4));
    
    
end

% pre-allocating arrays
C_low_tension = zeros(5,2);
C_mid_tension = zeros(5,2);
C_high_tension = zeros(5,2);

C_low_tension(1,1) = mean(C_bins(1,1:143));
C_low_tension(1,2) = mean(C_bins(2,1:143));

C_mid_tension(1,1) = mean(C_bins(1,144:177));
C_mid_tension(1,2) = mean(C_bins(2,144:177));

C_high_tension(1,1) = mean(C_bins(1,178:end));
C_high_tension(1,2) = mean(C_bins(2,178:end));

% writing saved initial concentrations in tension areas in .txt output files
csvwrite(['C_low_tension_' num2str(time) '.txt'],C_low_tension);
csvwrite(['C_mid_tension_' num2str(time) '.txt'],C_mid_tension);
csvwrite(['C_high_tension_' num2str(time) '.txt'],C_high_tension);

C_time_array = [0;28800;86400;259200;604800];

correction = zeros(1,200);
%--------------------------------------------------------------------------
while time < t_end

time = time + dt;  
for iii = 1:200
    
    % particle number in current bin iii
    ip = particle_array(:,2) == iii;
    
    % calculation of correction term
    if iii <= 199
         correction(iii) = ((((D_bins(iii)) - (D_bins(iii+1)))) / ((bin_edges(iii) - bin_edges(iii+1)) / 10^6)) * dt;      
    else      
         correction(iii) = ((((D_bins(iii-1)) - (D_bins(iii)))) / ((bin_edges(iii-1) - bin_edges(iii)) / 10^6)) * dt;
    end
    
    % calculation of displacement of these particles along LD
    delta_d = ((randn(sum(ip),1) * (sqrt(2 .* ((D_bins(iii))) .* dt))) - correction(iii))  * 10^6;
    
    % new position of particles
    particle_array(ip,1) = particle_array(ip,1) - delta_d;
    
    
end
%--------------------------------------------------------------------------

% IMPORTANT: BOUNDARIES. particles at the borders, which due to their
% random displacement left and right flying out of the domain,
% must be held in the domain AND ,according to their remainder
% displacement step, must be relocated back in the opposite direction to stay in domain 
while ~isempty(particle_array(particle_array(:,1) > bin_edges(1),1)) || ~isempty(particle_array(particle_array(:,1) < bin_edges(end),1))
particle_array(particle_array(:,1) > bin_edges(1),1) = bin_edges(1) - ((particle_array(particle_array(:,1) > bin_edges(1),1) - bin_edges(1)));
particle_array(particle_array(:,1) < bin_edges(end),1) = bin_edges(end) - (particle_array(particle_array(:,1) < bin_edges(end),1) - bin_edges(end));
end     
%--------------------------------------------------------------------------
% After being moved, particles get new bin numbers when they come into the range of a new bin,
% but keep their concentration. It is also possible that
% skip particle bins.
    new_bins_particles = (length(N_bins) + 1) - (discretize(particle_array(:,1),flip(bin_edges)));
    particle_array(:,2) = new_bins_particles;
    
%--------------------------------------------------------------------------
% calculation of concentrations in bins after displacement
for iii = 1:200
    
    ip = (particle_array(:,2) == iii);
    
    % mean conc. H in bins
    C_bins(1,iii) = mean(particle_array(ip,3));
    
    % mean conc. O in bins
    C_bins(2,iii) = mean(particle_array(ip,4));
    
    
end
%--------------------------------------------------------------------------
% saving of mean concentrations in three tension areas at certain points
% in time
if time == 28800 || time == 86400 || time == 259200 || time == 604800
   
    C_low_tension(C_time_array==time,1) = mean(C_bins(1,1:143));
    C_low_tension(C_time_array==time,2) = mean(C_bins(2,1:143));

    C_mid_tension(C_time_array==time,1) = mean(C_bins(1,144:177));
    C_mid_tension(C_time_array==time,2) = mean(C_bins(2,144:177));

    C_high_tension(C_time_array==time,1) = mean(C_bins(1,178:end));
    C_high_tension(C_time_array==time,2) = mean(C_bins(2,178:end));
    
    % writing saved concentrations in tension areas in .txt output files
    csvwrite(['C_low_tension_' num2str(time) '.txt'],C_low_tension);
    csvwrite(['C_mid_tension_' num2str(time) '.txt'],C_mid_tension);
    csvwrite(['C_high_tension_' num2str(time) '.txt'],C_high_tension);
end

end



