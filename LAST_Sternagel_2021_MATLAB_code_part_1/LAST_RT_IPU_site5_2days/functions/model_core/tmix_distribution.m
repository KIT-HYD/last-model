% Function to calculate time scale for mixing (tmix) of event and pre-event
% particles (with a possible distribution of mixing times)
function [t_mix_ratio, t_mix_mean1, t_mix_SD1, t_mix_mean2, t_mix_SD2] = tmix_distribution(D_table, D, dim, ks, dtc)

t_mix_ratio = 1.0; % to switch off tmix distribution set t_mix_ratio to 1.0. Then, t_mix_mean1 is the only mixing time.

icl = sum(D_table(1,:) < mean(D(1:dim-1)));
Dmix = prctile(D_table(1,icl+1:end),70);

t_mix_mean1 = ceil(((ks(1) .* dtc).^2) ./ Dmix);
t_mix_SD1 = 1;

t_mix_mean2 = 1;
t_mix_SD2 = 1;

end

