%% Function to read initial concentration profile of soil matrix to file

function [Cw_init]=readiniconc_state_interp(iniconc_file,z,dim)

%% Initialisation of parameters

init = dlmread(iniconc_file,' ');
[m_i , np] = size(init); % m_i number of observation values, np is number of columns (not necessary, just for information)
tiefe = zeros(m_i,1); % start depth
tiefe2 = zeros(m_i,1); % end depth
Cw_init = init(m_i,1)*ones(dim,1); % sets concentration in every depth to the value of concentration at lower boundary as initialisation
tiefe(1) = init(1,2);
tiefe2(1) = init(1,3);

%% Defines initial tracer concentration profile of soil matrix using interpolation

 for i=2:m_i
  tiefe(i) = init(i,2); % defines start depth for interpolation
  tiefe2(i) = init(i,3); % defines end depth for interpolation
  pos =find (z  <= -tiefe(i-1) & z  > -tiefe2(i-1)); % finds positions of z where the range between start and end depth is valid
  
      for j= 1:length(pos) % interpolation of concentration for every grid node of z
       Cw_init(pos(j))=init(i-1,1)+(init(i,1)-init(i-1,1))*(abs(z(pos(j)))-abs(tiefe(i-1)))/(abs(tiefe2(i-1))-abs(tiefe(i-1))); 
      end
  
 end

