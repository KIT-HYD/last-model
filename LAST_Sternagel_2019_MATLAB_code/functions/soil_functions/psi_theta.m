%% Function to calculate psi as function of theta

function [psi, c] = psi_theta(alph, dim, istyp, n_vg, stor, theta, thr, ths)

%% Initialisation of parameters

psi = zeros(dim,1);
c = zeros(dim,1);
y = zeros(dim,1);

%% Calculates initial psi and c state of soil matrix

for i = 1:dim
    m_vg = 1 - 1 ./ n_vg(istyp(i));
    th_star = (theta(i) - thr(istyp(i))) / (ths(istyp(i)) - thr(istyp(i)));    
    psi(i) = -1* ( (1 - th_star.^(1./m_vg))./ (th_star.^(1./m_vg)) ).^(1./n_vg(istyp(i))) ./alph(istyp(i));   
    
   
    if theta(i) <= 0.99 * ths(istyp(i))
       y(i) = -m_vg .* (1 ./ (1 + abs((psi(i)) .* alph(istyp(i))) .^ n_vg(istyp(i)))) .^ (m_vg+1) .* n_vg(istyp(i)) .* (abs(psi(i)) .* alph(istyp(i))) .^ (n_vg(istyp(i))-1) .* alph(istyp(i));
       c(i) = -(ths(istyp(i)) - thr(istyp(i))) .* y(i);
    elseif theta(i) > 0.99 * ths(istyp(i))
           c(i) = stor(istyp(i));
    end

end
