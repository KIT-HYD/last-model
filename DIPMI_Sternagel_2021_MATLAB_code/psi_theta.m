%% Function to calculate psi as function of theta

function [bins_psi] = psi_theta(alph, n_vg, theta, thr, ths)

%% psi in m Wassersäule

% psi = zeros(dim,1);

m = 1 - (1 ./ n_vg);
S = (theta - thr) / (ths - thr);    
% bins_psi = -1 * ((((1 - th_star.^(1./m_vg))./ (th_star.^(1./m_vg))).^(1./n_vg)) ./ alph);  
bins_psi = -(((S.^(-1/m)) - 1).^(1/n_vg)) / alph;

bins_psi(bins_psi == -Inf) = bins_psi(2);

% umrechnung_m_HPa = 0.00980638;
% 
% bins_psi = bins_psi * umrechnung_m_HPa;
