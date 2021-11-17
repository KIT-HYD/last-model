% Function to calculate reactive transport with adsorption and degradation
% in matrix
function [Cw_ret, Ct, Ct_deg, Ct0, Cm, position_z] = func_reactive_transport(dim, position_z, z, solubility, m, rho, Ct,...
                                                                        Ct_deg, dz, theta, KF, beta, koeff, dtc)
                                                            
Ct0 = zeros(dim,1);
Cm = zeros(dim,1);
Cw_ret = zeros(dim,1);
particles = zeros(dim,1);

for i = 1:dim-1 
   ipart = (position_z(:,1) <= z(i) & position_z(:,1) > z(i+1));
   particles(i) = sum(ipart); % determining number of particles in grid element
   
   % mixing of reactive solute masses among particles
   Cw_ret(i) = sum(position_z(ipart,4)) / particles(i);
   position_z(ipart,4) = Cw_ret(i) * ones(sum(ipart),1);
        
   if Cw_ret(i) > solubility % check for solubility
      Ct0(i) = sum((position_z(ipart,4) .* (m / rho)) - (solubility .* (m / rho))); % residual, insoluble amount of solute is stored in own array. this amount of solute is assmued to be precipitated in the actual grid element.
      Ct(i) = Ct(i) + Ct0(i);
      Ct_deg(i) = Ct_deg(i) + Ct0(i);
      position_z(ipart,4) = solubility * ones(sum(ipart),1);
      Cw_ret(i) = sum(position_z(ipart,4)) / particles(i);
   end

   %retardation/sorption
   %case1: normal, adsorption of solutes from water phase to adsorbing phase
    if (Ct_deg(i) / dz(i)) < (Cw_ret(i) * theta(i))           
      Cw_ret_particle = real((KF(i)/1000) .* ((position_z(ipart,4)).^beta)); % calculation of sorption, first part of Eq.6         
      m_ret_particle = sum(Cw_ret_particle .* (m / rho)); % conversion from conc. to masses, second part of Eq. 6
      Ct(i) = Ct(i) + m_ret_particle; % retarded/sorbed mass stored as solute mass (kg) in adsorbed phase (Ct = array containing adsorbed solutes masses, which are not degraded, only sorbed)
      Ct_deg(i) = Ct_deg(i) + m_ret_particle; % retarded/sorbed mass stored as solute mass (kg) in adsorbed phase (Ct_deg = array containing adsorbed solutes masses, which are then degraded, sorption + degradation)
       
      % check for minus conc. of particles
      if position_z(ipart,4) - Cw_ret_particle < 0
         position_z(ipart,4) = 0;
      else
         position_z(ipart,4) = position_z(ipart,4) - Cw_ret_particle;
      end

      Cw_ret(i) = sum(position_z(ipart,4)) / particles(i); % calculating retarded conc. profile in water phase

   %case2: reverse, desorption of solutes from adsorbing phase to water
   %phase, only when concentration in adsorbing phase is higher and
   %solubility in water phase allows to gather more solutes 
    elseif ((Cw_ret(i) * theta(i)) < (Ct_deg(i) / dz(i))) && round(Cw_ret(i),4) < round(solubility,4)          
      m_des = real((KF(i)/1000) .* ((Ct_deg(i) / dz(i)).^beta)) * dz(i); % desorbing mass from adsorbing phase to water phase
      
      if ((((particles(i) * m) / rho) * Cw_ret(i)) + m_des) > (((particles(i) * m) / rho) * solubility) % checks if addition of desorbed mass exceeds solubility in water phase of matrix
            m_des = m_des - (((((particles(i) * m) / rho) * Cw_ret(i)) + m_des) - (((particles(i) * m) / rho) * solubility));
      end
      
      Cw_des = (m_des / particles(i)) / (m / rho); % conversion from masses to conc.
      position_z(ipart,4) = position_z(ipart,4) + Cw_des; % adding desorbed conc. or mass to water phase and particles
      Cw_ret(i) = sum(position_z(ipart,4)) / particles(i); % updating of retarded conc. profile
      Ct(i) = Ct(i) - m_des; % updating mass stored in adsorbing phase
      Ct_deg(i) = Ct_deg(i) - m_des; % updating mass stored in adsorbing phase

    end

   %case3: balance, conc. of water phase and adsorbing phase are equal,
   % no sorption at all

   %degradation (only in adsorbed phase)
   Ct_deg(i) = real(Ct_deg(i) * (1 - koeff(i) * dtc/86400));
   Cm(i) = Ct(i) - Ct_deg(i); 

end



end

