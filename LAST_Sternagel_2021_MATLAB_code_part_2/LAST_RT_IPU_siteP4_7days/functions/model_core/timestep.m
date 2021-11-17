% Function to set adaptive, dynamic time step
function dtc = timestep(time, t_boundary, m_surface, m, n_mak, qb_u, pfd_z, pfd_dim, pfd_v, v, D, dz, prec_time, dtc)

if (m_surface > m && n_mak > 0) || (qb_u ~= 0 && n_mak > 0)
    dtc = ceil(abs(pfd_z(pfd_dim))/pfd_v); % when active infiltration into macropores the time step is set to small values
else
   % larger time step during active infiltration without macropores or generally inactive infiltration ("3.85" is not a random factor but an average value for a set of uniformly distributed numbers)
   while ((max(v) + (0.25 * max(diff(D)) / dz(1)) * dtc) + (3.85 .* (2 * max(D) * dtc) .^ 0.5)) / 0.1 < (1 * 0.99)
        dtc = dtc + 1;
   end
   % double while loop ensures that the Courant criterion is right beneath 1
   while ((max(v) + (0.25 * max(diff(D)) / dz(1)) * dtc) + (3.85 .* (2 * max(D) * dtc) .^ 0.5)) / 0.1 >= (1 * 0.99)
        dtc = dtc - 1;
   end 
end

if dtc > prec_time(2) - prec_time(1)
   dtc = prec_time(2) - prec_time(1); 
end

if (dtc + time) > prec_time(find(t_boundary <= time, 1, 'last' )+1)
    dtc = prec_time(find(t_boundary <= time, 1, 'last' )+1) - time;
end

end

