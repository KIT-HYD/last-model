%% Function to initialise particle positions and concentrations
% takes the number of particles in each grid element (n(i)) and defines each of
% these particles within a grid element (between z(i+1) and z(i) in z
% direction) with a certain location based on normal distribution

function [position_z, c_particle] = init_particles_pos(z, dim, n, Cw_init)

%% Initialisation of position and concentration arrays
% 
% position_z = cell(1,dim-1);
% c_particle = cell(1,dim-1);
% 
% %% Calculates positions and concentrations of every particle
% 
% for i=1:(dim-1)
%     increment = (abs(z(i))+z(i+1))/n(i);
%     position_z{i} = (z(i):increment:(z(i+1)-increment)); % normal distribution of particles in each grid element
%     c_particle{i} = (Cw_init(i)*ones(n(i),1))'; % every particle in the matrix gets the concentration of its grid element
% end
% 
% position_z = round([position_z{:}],3)';
% c_particle = round([c_particle{:}]',3);


for i=1:(dim-1)
    s_z(i)=arrayfun(@(x,y,p) x+(y-x)*rand(1,p),z(i+1),z(i),n(i),'UniformOutput', false);
end
position_z=[s_z{:}]'; % initial z position of the particles
c_particle = 0 .* ones(length(position_z),1);

                  