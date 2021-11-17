%% Function to calculate the mass balances and plot the mass profile

function [total_deg, recovery_sim, deg_recovery_sim,recovery_obs] = plot_mass(pfd_Ct_deg, pfd_Ct, Ct, Ct_deg, Cw_ret, prec_int, conc, concfin, theta, dim, Cw, dz, theta_init, Cw_init, z, z_plot )

%% Calculates mass balance

mass_input=sum(prec_int.*conc*(prec_int(2) - prec_int(1)))*1000; % total input mass
mass_obs_final=sum(concfin(1:length(concfin)-1,1).*(concfin(2:length(concfin),3)-concfin(2:length(concfin),2))*1*1)*1000; % final observed mass  g
mass_sim_final=sum(theta(1:dim-1).*Cw(1:dim-1).*dz)*1000; % final simulated mass g 
ret_mass_sim_final = sum(theta(1:dim-1).*(Cw_ret(1:dim-1)).*dz)*1000; % final simulated mass with retardation g
%total_ret = mass_sim_final - ret_mass_sim_final; % total retarded mass g
total_deg = sum(Ct - Ct_deg) * 1000 + sum(pfd_Ct - pfd_Ct_deg) * 1000; % total degraded mass g

recovery_sim=(mass_sim_final/mass_input)*100; % recovery rate of simulated masses
deg_recovery_sim = ((ret_mass_sim_final + sum(Ct_deg(1:dim-1)*1000))/mass_input)*100; % recovery rate of degraded masses
recovery_obs= (mass_obs_final/mass_input)*100; % recovery rate of observed masses

%% Mass plot
z_plot = z;
figure;

hold on;
% h1=plot((theta_init(1:dim-1).*Cw_init(1:dim-1).*dz*1*1)*1000,z_plot(1:dim-1)','k-','linewidth',2,'markersize',2);
h2=plot((theta(1:dim-1).*Cw(1:dim-1).*dz*1*1)*1000,z_plot(1:dim-1)','r-','linewidth',2,'markersize',2);
h3=plot((concfin(:,1).*0.1*1*1)*1000,-concfin(:,2),'b-','linewidth',2,'markersize',2);
h4=plot(((theta(1:dim-1).*abs(Cw_ret(1:dim-1)).*dz*1*1)*1000) + (Ct(1:dim-1)*1000),z_plot(1:dim-1)','y-','linewidth',2,'markersize',2);
% h5=plot((theta(1:dim-1).*abs(Cw_ret(1:dim-1)).*dz*1*1)*1000,z_plot(1:dim-1)','g-','linewidth',2,'markersize',2);
h6=plot(((theta(1:dim-1).*abs(Cw_ret(1:dim-1)).*dz*1*1)*1000) + (Ct_deg(1:dim-1)*1000),z_plot(1:dim-1)','c-','linewidth',2,'markersize',2);
% h7=plot(Ct_deg(1:dim-1)*1000,z_plot(1:dim-1)','k--','linewidth',2,'markersize',2);
hold off;
title(['Solute Mass Profile at End of Simulation (Input=' num2str(mass_input) 'g)'],'fontsize',14);
xlabel('Mass [g]','fontsize',14);
ylabel('z [m]','fontsize',14);
legend([h2 h3 h4 h6],'conservative' ,'Final Mass Profile (Observed)','only sorption','sorption+degradation','Location','southeast');
ylim([-1.0 0 ]);
xlim([0 0.5]);
set(gca,'fontsize',14,'linewidth',2,'XMinorTick','on','YMinorTick','on');


end

