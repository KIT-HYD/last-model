%% Function to calculate the mass balances and plot the mass profile

function plot_mass(prec_int, conc, concfin, theta, dim, Cw, dz, theta_init, Cw_init, z, z_plot )

%% Calculates mass balance

mass_input=sum(prec_int.*conc*600)*1000; % total input mass
mass_obs_final=sum(concfin(:,1).*(concfin(:,3)-concfin(:,2))*1.4*1.4)*1000; % final observed mass
mass_sim_final=sum(theta(1:dim-1).*Cw(1:dim-1).*dz)*1000; % final simulated mass

recovery_sim=(mass_sim_final/mass_input)*100; % recovery rate of simulated masses
recovery_obs= (mass_obs_final/mass_input)*100; % recovery rate of observed masses

%% Mass plot

figure;

hold on;
h1=plot((theta_init(1:dim-1).*Cw_init(1:dim-1).*dz*1.4*1.4)*1000,z_plot(1:dim-1)','k-','linewidth',2,'markersize',2);
h2=plot((theta(1:dim-1).*Cw(1:dim-1).*dz*1.4*1.4)*1000,z_plot(1:dim-1)','r-','linewidth',2,'markersize',2);
h3=plot((concfin(:,1).*(concfin(:,3)-concfin(:,2))*1.4*1.4)*1000,-concfin(:,3),'b-','linewidth',2,'markersize',2);
hold off;
title(['Solute Mass Profile at End of Simulation (Input=' num2str(mass_input) 'g)'],'fontsize',14);
xlabel('Mass [g]','fontsize',14);
ylabel('z [m]','fontsize',14);
legend([h1 h2 h3],'Initial Mass Profile',['Final Mass Profile (Simulated), Mass =' num2str(mass_sim_final) 'g , Recovery =' num2str(recovery_sim) '%'] ,['Final Mass Profile (Observed), Mass =' num2str(mass_obs_final) 'g , Recovery =' num2str(recovery_obs) '%'],'Location','southeast');
ylim([-1.0 0 ]);
xlim([0 mass_input]);
set(gca,'fontsize',14,'linewidth',2,'XMinorTick','on','YMinorTick','on');


end

