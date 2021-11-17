%% Function to read initial theta or psi state from file

function [theta_init, psi_init] = readini_state_interp(inistate_file,i_init,z,dim)

%% Initialisation of parameters

init = dlmread(inistate_file,' ');
[m_i , np] = size(init); % m_i number of observation values, np is number of columns (not necessary, just for information)
tiefe = zeros(m_i,1); % start depth
tiefe2 = zeros(m_i,1); % end depth
theta_init = init(m_i,1) * ones(dim,1); % sets soil moisture in every depth to the value of moisture at lower boundary as initialisation
psi_init = zeros(dim,1); % sets matrix potential in every depth to 0 as initialisation
tiefe(1) = init(1,2);
tiefe2(1) = init(1,3);

%% Defines initial theta or psi conditions of soil matrix dependent on input type using interpolation

if i_init == 1 % theta as input type
  
     for i=2:m_i
      tiefe(i) = init(i,2); % defines start depth for interpolation
      tiefe2(i) = init(i,3); % defines end depth for interpolation
      pos = find (z  <= -tiefe(i-1) & z  > -tiefe2(i-1)); % finds positions of z where the range between start and end depth is valid

          for j=1:length(pos) % interpolation of soil moisture for every grid node of z
           theta_init(pos(j)) = init(i-1,1) + (init(i,1)-init(i-1,1)) * (abs(z(pos(j))) - abs(tiefe(i-1))) / (abs(tiefe2(i-1)) - abs(tiefe(i-1))); 
          end

     end
 
elseif i_init == 2 % psi as input type
    
     for i=1:m_i
      psi_s    = init(i,1);
      tiefe(i) = init(i,2);
      tiefe2(i) = init(i,3);
      pos = find (z  <= -tiefe(i) & z  > -tiefe2(i));
      
          for j= 1:length(pos)
           psi_init(pos(j)) = psi_s;
          end
      
     end
end

