%% Function to calculate theta as function of psi

function [theta, c]=theta_psi(alph,dim, istyp, n_vg, psi, stor, thr, ths)

%% Initialisation of parameters

v2struct(soil_ws);
theta=zeros(dim,1);
c=zeros(dim,1);
y=zeros(dim,1);

%% Calculates initial theta and c state of soil matrix

for i=1:dim
    
       if psi(i)>=0
         psi(i)=-0.001;
       end
   
    m_vg = 1 - 1./n_vg(istyp(i));
    th_star = (1./(1+(abs(psi(i)).*alph(istyp(i))).^n_vg(istyp(i)))).^m_vg;
    theta(i) = thr(istyp(i))+th_star*(ths(istyp(i))-thr(istyp(i)));
   
       % check constraints
       if theta(i) > ths(istyp(i))
          theta(i)= ths(istyp(i));
       elseif theta(i) < thr(istyp(i))
              theta(i) = thr(istyp(i))*1.1;
       end 
       
       if theta(i) <= 0.99* ths(istyp(i))
          y(i)= -m_vg.*(1./(1+abs((psi(i)).*alph(istyp(i))).^n_vg(istyp(i)))).^(m_vg+1) .*n_vg(istyp(i)).*(abs(psi(i)).*alph(istyp(i))).^(n_vg(istyp(i))-1).*alph(istyp(i));
          c(i)=-(ths(istyp(i))-thr(istyp(i))).*y(i);
       elseif theta(i) > 0.99 *ths(istyp(i))
              c(i)=stor(istyp(i));
       end
end
