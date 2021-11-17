%% Function for the calculation of the pre-event particles displacement

function [Cw, mtx_avg_age, position_znew, theta] = ... % output arguments
            particle_track_new_mob... % function
            (D, dim, dt, D_table, dz, istyp, K_table, m, mob_fak, position_z, ths, v, z) % input arguments

%% Initialisation of parameters and arrays

position_znew = position_z; % position of every particle in soil matrix
particles = zeros(dim,1); % amount of particles in each grid element of soil matrix
theta = zeros(dim,1); % soil moisture profile of soil matrix
Cw = zeros(dim,1); % concentration profile of soil matrix
mtx_avg_age = zeros(dim,1); % age profile of soil matrix  
M_layer = zeros(dim,1); % total water mass in each layer
ip = zeros(dim,1); % storage for particle amounts in each grid element

ip_0 = position_z(:,1) > z(1); % correction for unphysical particle losses at top of the soil matrix
position_z(ip_0) = z(1) - position_z(ip_0); % ensures that every particle has a negative position value within the soil matrix

%% Particle displacement

for j = 1:dim-1
   icl = sum(D_table(istyp(j)) < D(j)); % finds the last bin where the assumed diffusivity of the lookup table (D_table) is still smaller than the real diffusivity (D). 'icl' is then the number of the bin with the maximum possible diffusivity 
   ip(j+1) = sum(position_z(:,1) <= z(j) & position_z(:,1) > z(j+1)); % finds the amount of all particles within the actual grid element (between z(j) and z(j+1))
   
   if ip(j+1) > 0
     d_part = floor(ip(j+1)/(icl+1)); % allocates the total amount of particles within a grid element to the possible amount of diffussivity bins. 'd_part' is then the total amount of particles within a diffussivity bin
     mobile = floor(mob_fak*ip(j+1)); % fraction of particles which contribute to the particle movement, thus the "mobile" fraction
     
     % defining of scaling factors
     mo_fak = zeros(ip(j+1),1); % initialisation of factor for scaling diffussive motion of particles in different bins
     k_fak = zeros(ip(j+1),1);
     for ii=1:icl % establishes a factor for scaling the diffusivity of the particles stored in the different diffussivity bins
         mo_fak((ii-1)*d_part+1:ii*d_part) = D_table(istyp(j),ii) / D(j) * ones(d_part,1);  % 'mo_fak' is thereby just the relation of the diffusivity of the single bins to the real diffussivity in the actual grid element. This factor is then applied to every particle within the bins
         k_fak((ii-1)*d_part+1:ii*d_part) = K_table(istyp(j),ii) / v(j) * ones(d_part,1); % 'k_fak' is thereby just the relation of the hydraulic conductivity of the single bins to the real hydraulic conductivity in the actual grid element. This factor is then applied to every particle within the bins
     end
     mo_fak(icl*d_part+1:length(mo_fak)) = ones(length(mo_fak)-icl*d_part,1); % these particles not considered by the floor function just get the highest scaling factor of one
     k_fak(icl*d_part+1:length(k_fak)) = ones(length(k_fak)-icl*d_part,1); % these particles not considered by the floor function just get the highest scaling factor of one

     % interpolates v and D for all particles within the grid element
     [v_part, D_part] = interpolation(v(j), v(j+1), D(j), D(j+1), z(j), z(j+1), position_z((sum(ip)-ip(j+1))+1)); 

     % calculates the random walk term
     correction = 0.25 * (D(j)-D(j+1)) .* mo_fak(ip(j+1)-mobile+1:ip(j+1)) / dz(1); % correction term for random walk
     randomwalk = 1 .* 1 .* ((randn(mobile,1))) .* mo_fak(ip(j+1)-mobile+1:ip(j+1)) .* ((2*D_part*dt) .^ 0.5); % random walk term
     %randomwalk = 1.0e-04 * ones(length(correction),1);
     position_znew((sum(ip)-mobile)+1:sum(ip),1) = real(position_z((sum(ip)-mobile)+1:sum(ip),1)... % calculates the new positions of the particles in each grid element
       - ((1*v_part*k_fak(ip(j+1)-mobile+1:ip(j+1))+correction) * dt + randomwalk)); 

   end
end

% postprocessing for realistic particle positions
ip_0 = position_znew(:,1) >= z(1); % correction for unphysical particle losses at top of the soil matrix due to the displacement process
position_znew(ip_0) = z(1) - position_znew(ip_0);

ip_low = position_znew(:,1) <= z(dim); % correction for unphysical particle losses at bottom of the soil matrix due to the displacement process
position_znew(ip_low) = z(dim-1); 


%% Calculation of water and solute masses in each grid element and generating of new soil moisture and concentration profile


for i=1:dim-1
    ipart = (position_znew(:,1) <= z(i) & position_znew(:,1) > z(i+1));
    particles(i) = sum(ipart); % amount of particles in grid element
    
    % check for oversaturation
    if particles(i) > round(ths(istyp(i))*dz(i)*1*1000/m) %if the number of particles is higher than the saturated threshold the number is set to this threshold
       particles(i)= round(ths(istyp(i))*dz(i)*1*1000/m);
       
       ipart2 = find(position_znew(:,1) <= z(i) & position_znew(:,1) > z(i+1));
       for ii=round(ths(istyp(i))*dz(i)*1*1000/m)+1:length(ipart2) 
           position_znew(ipart2(ii),1)= position_znew(ipart2(ii),1) - dz(i);
       end
    end
    
    % calculation of new profiles
    M_layer(i) = particles(i) * m; % total mass in each grid element
    theta(i) = M_layer(i) / (dz(i)*1*1000); % new soil moisture in each grid element
    
    mtx_avg_age(i) = sum(position_znew(ipart,2)) / particles(i); % new average age in each grid element
    if isnan(mtx_avg_age(i))
       mtx_avg_age(i) = 0;
    end 
    
    Cw(i)= sum(position_znew(ipart,3)) / particles(i); % new concentration in each grid element
    position_znew(ipart,3) = Cw(i) * ones(sum(ipart),1); % new concentration of every single particle in each grid element
                
end

end