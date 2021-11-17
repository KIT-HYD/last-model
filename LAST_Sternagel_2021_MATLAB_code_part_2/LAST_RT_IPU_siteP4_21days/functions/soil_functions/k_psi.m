%% Function to calculate soil hydraulic conductivity from psi

function k = k_psi(alph, dim, istyp, ks, l_vg, n_vg, psi)

%% Initialisation of parameters

k = zeros(dim,1);

%% Calculates initial hydraulic conductivities of soil matrix

  for i=1:dim
    m_vg = 1 - 1./n_vg(istyp(i));
    v = 1 + (alph(istyp(i)) .* abs(psi(i))).^n_vg(istyp(i));
    k(i) = ks(istyp(i)) .* (1 - (1 - v.^(-1)).^m_vg).^2 ./ (v.^(m_vg.*l_vg));
   end 
