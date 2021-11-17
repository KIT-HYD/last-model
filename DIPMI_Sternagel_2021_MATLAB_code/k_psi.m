%% Function to calculate soil hydraulic conductivity from psi

function k = k_psi(alph, dim, ks, l_vg, n_vg, psi)

%% Initialisation of parameters

k = zeros(dim,1);

%% Calculates initial hydraulic conductivities of soil matrix

  for i=1:dim
    m_vg = 1 - 1./n_vg;
    v = (1 + (alph .* abs(psi(i))).^n_vg);
    k(i) = ks .* (1 - (1 - v.^(-1)).^m_vg).^2 ./ (v.^(m_vg.*l_vg));
   end 