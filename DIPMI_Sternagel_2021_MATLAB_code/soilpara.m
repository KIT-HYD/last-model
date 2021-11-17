%% Function to read soil parameter file
function [alph, ks, l_vg, n_vg, stor, thr, ths] = soilpara(soilpara_file)

%% Initialisation of parameters

para = dlmread(soilpara_file,' ');
[m_s , np] = size(para); % m_s number of soil types, np is number of soil parameters (not necessary, just for information)
thr   =  zeros(m_s,1); % residual soil moisture
ths   =  zeros(m_s,1); % soil moisture at saturation
alph =  zeros(m_s,1); % bubbling pressure
n_vg  =  zeros(m_s,1); % n (width of pore size distribution)
ks    =  zeros(m_s,1); % saturated hydraulic conductivity
stor  =  zeros(m_s,1); % specific storage capacity
l_vg  = 0.5; 
% istyp = m_s * ones(dim,1); % pointer on soil types (number of different soil types in the domain (1= just one soil typ all over the whole depth))
tiefe = zeros(m_s,1); % start depth
tiefe2 = zeros(m_s,1); % end depth

%% Defines parameters for all soil types

for i=1:m_s
 ks(i)    = para(i,1);
 ths(i)   = para(i,2);
 thr(i)   = para(i,3); 
 alph(i) = para(i,4); 
 n_vg(i)  = para(i,5);
 stor(i)  = para(i,6);
 tiefe(i) = para(i,7);
 tiefe2(i) = para(i,8);
 
% % looks for depths where actual soil type is valid
% pos = find (z  <= -tiefe(i) & z  > -tiefe2(i)); % for just one soil type => every position of z between start depth (0m) and end depth (i.e. -1,5m)
% 
%     for j= 1:length(pos)
%     istyp(pos(j)) = i;
%     end

end
