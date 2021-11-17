% Function to calculate reactive transport with adsorption and degradation
% in pfd
function [pfd_Cw_ret, pfd_Ct, pfd_Ct_deg, pfd_Ct0, pfd_Cm, pfd_position_z] = func_pfd_reactive_transport(pfd_dim, pfd_position_z, pfd_z, solubility, pfd_m, rho, pfd_Ct,...
                                                                        pfd_Ct_deg, pfd_dz, pfd_theta, pfd_KF, beta, pfd_koeff, dtc)
                                                            
pfd_Ct0 = zeros(pfd_dim,1);
pfd_Cm = zeros(pfd_dim,1);
pfd_Cw_ret = zeros(pfd_dim,1);
particles = zeros(pfd_dim,1);

for i=1:pfd_dim-1 
  if pfd_theta(i) > 0  
   ipart = (pfd_position_z(:,1) <= pfd_z(i) & pfd_position_z(:,1) > pfd_z(i+1));
   particles(i) = sum(ipart);
        
   pfd_Cw_ret(i) = sum(pfd_position_z(ipart,5)) / particles(i);
   pfd_position_z(ipart,5) = pfd_Cw_ret(i) * ones(sum(ipart),1);
        
%    if pfd_Cw_ret(i) > solubility % check for solubility
%       pfd_Ct0(i) = sum((pfd_position_z(ipart,5) .* (pfd_m / rho)) - (solubility .* (pfd_m / rho))); % residual, insoluble amount of solute is stored in own array. this amount of solute is assmued to be precipitated in the actual grid element.
%       pfd_Ct(i) = pfd_Ct(i) + pfd_Ct0(i);
%       pfd_Ct_deg(i) = pfd_Ct_deg(i) + pfd_Ct0(i);
%       pfd_position_z(ipart,5) = solubility * ones(sum(ipart),1);
%       pfd_Cw_ret(i) = sum(pfd_position_z(ipart,5)) / particles(i);
%    end

   %retardation
   %case1: normal, adsorption of solutes from water phase to adsorbing phase
    if (pfd_Ct_deg(i) / pfd_dz(i)) < (((particles(i) * pfd_m) / rho) * pfd_Cw_ret(i)) / pfd_dz(i)          
      Cw_ret_particle = real((pfd_KF(i)/1000) .* ((pfd_position_z(ipart,5)).^beta));          
      m_ret_particle = sum(Cw_ret_particle .* (pfd_m / rho)); % conversion from conc. to masses
      pfd_Ct(i) = pfd_Ct(i) + m_ret_particle; % retarded conc. or mass stored as solute mass kg in adsorbed phase
      pfd_Ct_deg(i) = pfd_Ct_deg(i) + m_ret_particle; % same for degradation array

      if pfd_position_z(ipart,5) - Cw_ret_particle < 0
         pfd_position_z(ipart,5) = 0;
      else
         pfd_position_z(ipart,5) = pfd_position_z(ipart,5) - Cw_ret_particle;
      end

      pfd_Cw_ret(i) = sum(pfd_position_z(ipart,5)) / particles(i); % calculating of retarded conc. profile

   %case2: reverse, desorption of solutes from adsorbing phase to water phase
    elseif ((((particles(i) * pfd_m) / rho) * pfd_Cw_ret(i)) / pfd_dz(i) < (pfd_Ct_deg(i) / pfd_dz(i))) % && round(pfd_Cw_ret(i),4) < round(solubility,4)          
      m_des = real((pfd_KF(i)/1000) .* ((pfd_Ct_deg(i) / pfd_dz(i)).^beta)) * pfd_dz(i); % desorbing mass from adsorbing phase to water phase
      
%       if ((((particles(i) * pfd_m) / rho) * pfd_Cw_ret(i)) + m_des) > (((particles(i) * pfd_m) / rho) * solubility)
%             m_des = m_des - (((((particles(i) * pfd_m) / rho) * pfd_Cw_ret(i)) + m_des) - (((particles(i) * pfd_m) / rho) * solubility));
%       end
      
      Cw_des = (m_des / particles(i)) / (pfd_m / rho); % conversion from masses to conc.
      pfd_position_z(ipart,5) = pfd_position_z(ipart,5) + Cw_des; % adding desorbed conc. or mass to water phase and particles
      pfd_Cw_ret(i) = sum(pfd_position_z(ipart,5)) / particles(i); % updating of retarded conc. profile
      pfd_Ct(i) = pfd_Ct(i) - m_des; % updating mass stored in adsorbing phase
      pfd_Ct_deg(i) = pfd_Ct_deg(i) - m_des; % same for degradation array

    end

   %case3: balance, conc. of water phase and adsorbing phase are equal,
   % no sorption at all
  end
  
   %degradation (only in adsorbed phase)
   pfd_Ct_deg(i) = real(pfd_Ct_deg(i) * (1 - pfd_koeff(i) * dtc/86400));
   pfd_Cm(i) = pfd_Ct(i) - pfd_Ct_deg(i); 

  
end



end

