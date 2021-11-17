%% Function to calculate psi as function of theta

function [psi, c] = c_theta(alph, dim, n_vg, stor, theta, thr, ths)

%% Initialisation of parameters

psi = zeros(dim,1);
c = zeros(dim,1);
y = zeros(dim,1);

%% Calculates initial psi and c state of soil matrix

for i = 1:dim
    m_vg = 1 - 1 ./ n_vg;
    th_star = (theta(i) - thr) / (ths - thr);    
%     psi(i) = real(-1* ( (1 - th_star.^(1./m_vg))./ (th_star.^(1./m_vg)) ).^(1./n_vg) ./alph);   
    psi(i) = real(-(((th_star.^(-1/m_vg)) - 1).^(1/n_vg)) / alph);
   
    if theta(i) <= 0.99 * ths
       y(i) = -m_vg .* (1 ./ (1 + abs((psi(i)) .* alph) .^ n_vg)) .^ (m_vg+1) .* n_vg .* (abs(psi(i)) .* alph) .^ (n_vg-1) .* alph;
       c(i) = -(ths - thr) .* y(i);
    elseif theta(i) > 0.99 * ths
           c(i) = stor;
    end

end