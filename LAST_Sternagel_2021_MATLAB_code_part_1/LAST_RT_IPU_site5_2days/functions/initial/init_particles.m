%% Function to initialise particle mass and number in each grid element

function [n, m]=init_particles(theta, dz, dim, N)

%% Initialisation of parameters

M = ones(dim-1,1);
n = ones(dim-1,1);

%% Calculates particles mass based on soil water content in grid elements
% distributes the particles over the depth according to the soil moisture

for i=1:(dim-1) 
    M(i) = theta(i)*dz(i)*1*1000; % mass stored in each layer (kg)
    n(i) = round((theta(i)/(sum(theta))*N)); % particle number in each layer 
end

m = sum(M) / N; % mass of one particle
theta_p = n*m / (dz(1)*1*1000); %control calculation of soil moisture in each layer
m = theta(1) / theta_p(1)*m; % error calculation of soil moisture to correct the value for the mass of one particle (kg)