%% Function for describing the mixing between event and pre-event particles within the first grid elements of the soil matrix

function [c_particle, Cw_event, position_z, position_z_event] = mixing_mtx(age_event, c_particle, Cw_event, position_z, position_z_event, time, t_mix)

ip_mix=find(age_event < time-t_mix); % finds the positions of all particles which are old enough to contribute to the mixing process at given mixing time and elapsed simulation time
ip_res=find(age_event >= time-t_mix); % finds the positions of the residual particles which are still to young for mixing
   
    if ~isempty(ip_mix)
        ipp = find(ip_mix < length(position_z_event)); % selects from all particles just these which are able to contribute to mixing
        
        % adds the mixing particles to the position and concentration
        % arrays of the pre-event particles and deletes them out of the
        % event arrays
        position_z = [position_z_event(ip_mix(ipp));position_z]; 
        c_particle = [Cw_event(ip_mix(ipp));c_particle]; 
        position_z_event = position_z_event(ip_res); 
        Cw_event = Cw_event(ip_res);
    end  

end

