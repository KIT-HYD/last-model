%% Function for describing the mixing between event and pre-event particles within the first grid elements of the soil matrix

function [position_z, position_z_event] = mixing_mtx(dim, z, reactive_transport,position_z, position_z_event, time)

ip_mix = find(position_z_event(:,2) < time - position_z_event(:,3)); % finds the positions of all particles which are old enough to contribute to the mixing process at given mixing time and elapsed simulation time
ip_res = position_z_event(:,2) >= time - position_z_event(:,3); % finds the positions of the residual particles which are still to young for mixing
   
    if ~isempty(ip_mix)
        ipp = find(ip_mix <= length(position_z_event(:,1))); % selects from all particles just these which are able to contribute to mixing
        
        % adds the mixing particles to the position and concentration
        % arrays of the pre-event particles and deletes them out of the
        % event arrays
        if reactive_transport == true
            pos_age_conc_mix = [position_z_event(ip_mix(ipp),1), position_z_event(ip_mix(ipp),2), position_z_event(ip_mix(ipp),5), position_z_event(ip_mix(ipp),4)];         
        else
            pos_age_conc_mix = [position_z_event(ip_mix(ipp),1), position_z_event(ip_mix(ipp),2), position_z_event(ip_mix(ipp),5)];
        end 
        
        position_z = [pos_age_conc_mix;position_z];
        position_z_event = position_z_event(ip_res,:);
        
        for i = 1:dim-1
            ip = find(position_z(1:length(pos_age_conc_mix),1) <= z(i) & position_z(1:length(pos_age_conc_mix),1) > z(i+1));
            position_z(ip,1) = -(z(i) + (z(i) - z(i+1)) .* rand(length(ip),1));                     
        end
              
    end  

end

