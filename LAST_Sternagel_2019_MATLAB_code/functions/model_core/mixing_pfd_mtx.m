%% Function for describing the mixing between the pfd and the soil matrix based on the macropore depth distribution

function [c_particle, pfd_age, pfd_Cw_event, pfd_position_z, position_z] = ... % output arguments
            mixing_pfd_mtx... % function
                (c_particle, dtc, k, ks, m, mak_sml, mak_mid, mtx_particles_mix0, n_mak, particle_contact_grid,... % input arguments
                pfd_age, pfd_Cw_event, pfd_dim, pfd_dz, pfd_m, pfd_n, pfd_particles, pfd_particle_output,...
                pfd_position_z, pfd_r, pfd_theta, pfd_z, position_z, psi, rate_big, rate_mid, rate_sml, rho, theta, ths, z)


    for yy = pfd_dim-1:-1:1  

    %% Finds the respective depths of the soil matrix corresponding to the depths of the three macropore sizes

    indx_big = find(pfd_z(yy) > round(z,2)); % for big macropore

        if (yy-((pfd_dim-1)-find(round(pfd_z,2)==mak_mid))-1) > 0 % for medium macropores
           indx_mid = find(pfd_z(yy-((pfd_dim-1)-find(round(pfd_z,2)==mak_mid))-1) > round(z,2));
        else
           indx_mid = NaN;
        end    

        if (yy-((pfd_dim-1)-find(round(pfd_z,2)==mak_sml))-1) > 0 % for small macropores
            indx_sml = find(pfd_z(yy-((pfd_dim-1)-find(round(pfd_z,2)==mak_sml))-1) > round(z,2));
        else
            indx_sml = NaN;
        end

    %% Calculates mixing mass and converts it into integer number of particles entering the soil matrix

    % assumption: diffusive mixinng only from saturated grid elements of pfd
    % into soil matrix, no back diffusion from matrix into pfd
        if pfd_theta(yy) == 1 && (theta(indx_big(1))+theta(indx_big(1)-1))/2 < ths 

        % calculates the flux density and the mass of diffusive mixing based on psi
        % and converts the mass to integer number of particles leaving the pfd 
        q_mix_psi = abs  ((2*ks*((k(indx_big(1))+k(indx_big(1)-1))/2))/(ks+((k(indx_big(1))+k(indx_big(1)-1))/2))*(((psi(indx_big(1))+psi(indx_big(1)-1))/2)/(pfd_r*2))); %diffussive flow and mass calculated with psi_init
        m_mix_psi = (q_mix_psi*(2*pi*pfd_r*pfd_dz(yy))) * rho * dtc * n_mak; %calculation of total mass leaving pfd in current depth
        pfd_particles_mix = floor(m_mix_psi/pfd_m); %amount particles leaving macropore

        % calculates the number of particles having contact to the lateral surface
        % of a pfd grid element and ensures that only the maximum possible number
        % of contact particles is leaving
        rate = (particle_contact_grid * n_mak)/(pfd_n(yy) * n_mak); 
        pfd_particles_contact = floor(pfd_particles(yy) * rate); 

            if pfd_particles_mix > pfd_particles_contact % maximum possible contact particles
               pfd_particles_mix = pfd_particles_contact;
            end

            if(pfd_particles_mix+pfd_particle_output) > pfd_particles(yy) % maximum possible particles still present in pfd grid element after draingae
               pfd_particles_mix = pfd_particles - pfd_particle_output;
            end

        % converts leaving pfd particles into integer number of matrix particles
        % entering the soil matrix
        mtx_particles_mix = floor((pfd_particles_mix*pfd_m)/m); 

        pfd_particles_mix = round(pfd_particles_mix-((((pfd_particles_mix*pfd_m)/m)-mtx_particles_mix)*m)/pfd_m);  % residual amount of particles not considered by "floor" function are recalculated to the pfd
            if yy == pfd_dim-1
               pfd_particles_mix = round(pfd_particles_mix-((((pfd_particle_output*pfd_m)/m)-mtx_particles_mix0)*m)/pfd_m); 
            end

        %% Deletion of leaving particle amount out of the pfd arrays

        indx2 = find(pfd_position_z(:,1) < pfd_z(yy) & pfd_position_z(:,1) >= pfd_z(yy+1)); 

        pfd_position_z(indx2(end)-pfd_particles_mix:indx2(end),:) = [];
        pfd_age(((indx2(end)-(pfd_particles_mix):indx2(end)))) = [];
        pfd_conc_output2 = sum(pfd_Cw_event(indx2(end)-pfd_particles_mix:indx2(end)))/pfd_particles_mix; % solute concentration of leaving particles
        pfd_Cw_event(indx2(end)-pfd_particles_mix:indx2(end)) = [];

        %% Adds the amount of particle entering the soil matrix to the arrays in the respective depths

        % determines how the total number of entering particles is distributed to
        % the different soil matrix depths based on the actual grid elements of the
        % three macropore sizes
        pos_mid = [];
        pos_sml = [];

            if ~isnan(indx_sml) % valid when all macropores are still unsaturated
               pos_sml = round(z(indx_sml(1)) + (z(indx_sml(1)-1)-z(indx_sml(1))) .* rand((round(mtx_particles_mix*rate_sml)),1),3);
               pos_mid = round(z(indx_mid(1)) + (z(indx_mid(1)-1)-z(indx_mid(1))) .* rand((round(mtx_particles_mix*rate_mid)),1),3);
               pos_big = round(z(indx_big(1)) + (z(indx_big(1)-1)-z(indx_big(1))) .* rand((round(mtx_particles_mix*rate_big)),1),3);
            elseif ~isnan(indx_mid) % valid when just big and medium macropores are still unsaturated
               pos_big = z(indx_big(1)) + (z(indx_big(1)-1)-z(indx_big(1))) .* rand((floor(mtx_particles_mix*((rate_sml/2)+rate_big))),1);
               pos_mid = z(indx_mid(1)) + (z(indx_mid(1)-1)-z(indx_mid(1))) .* rand((ceil(mtx_particles_mix*((rate_sml/2)+rate_mid))),1);
            else % valid when just big macropores are still unsaturated
               pos_big = z(indx_big(1)) + (z(indx_big(1)-1)-z(indx_big(1))) .* rand((round(mtx_particles_mix*(rate_big+rate_mid+rate_sml))),1);
            end  

        % adding the positions and concentrations of the new particle to the arrays
%         position_z = sort([pos_big;pos_mid;pos_sml;position_z],'descend');
        position_z = [pos_big;pos_mid;pos_sml;position_z];

        down2_big = sum(position_z <= z(indx_big(1)-1));
        c_particle = [c_particle(1:(length(position_z)-(down2_big)));pfd_conc_output2*ones(length(pos_big),1);c_particle((length(position_z)-(down2_big)+1):length(c_particle))];

            if ~isempty(pos_mid)  
                down2_mid = sum(position_z <= z(indx_mid(1)-1));
                c_particle = [c_particle(1:(length(position_z)-(down2_mid)));pfd_conc_output2*ones(length(pos_mid),1);c_particle((length(position_z)-(down2_mid)+1):length(c_particle))]; 
            end

            if ~isempty(pos_sml)
               down2_sml = sum(position_z <= z(indx_sml(1)-1));
               c_particle = [c_particle(1:(length(position_z)-(down2_sml)));pfd_conc_output2*ones(length(pos_sml),1);c_particle((length(position_z)-(down2_sml)+1):length(c_particle))]; 
            end


        end

    end

end

