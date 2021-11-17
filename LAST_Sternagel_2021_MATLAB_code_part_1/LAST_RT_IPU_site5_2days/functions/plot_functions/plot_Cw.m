%% Function to plot the final concentration profile

function plot_Cw(dz, Ct_deg, Ct, Cw_ret, concfin, Cw, Cw_init, dim, theta, theta_init, z, z_plot)

figure;

hold on;
h1=plot(Cw_init(1:dim-1).*theta_init(1:dim-1),z_plot(1:dim-1)','k-','linewidth',2,'markersize',2);
h2=plot((Cw(1:dim-1).*theta(1:dim-1)),z_plot(1:dim-1)','r-','linewidth',2,'markersize',2);
h3=plot(concfin(:,1),-concfin(:,2),'b-','linewidth',2,'markersize',2);
h4=plot(Ct(1:dim-1)./dz,z_plot(1:dim-1)','y-','linewidth',2,'markersize',2);
h5=plot((Cw_ret(1:dim-1).*theta(1:dim-1)),z_plot(1:dim-1)','g-','linewidth',2,'markersize',2);
h6=plot((Cw_ret(1:dim-1).*theta(1:dim-1))+(Ct_deg(1:dim-1)./dz),z_plot(1:dim-1)','c-','linewidth',2,'markersize',2);
h7=plot(Ct_deg(1:dim-1)./dz,z_plot(1:dim-1)','k--','linewidth',2,'markersize',2);
hold off;
title('Solute Concentration Profile at End of Simulation','fontsize',14);
xlabel('Concentration in Soil Phase [kg/m³]','fontsize',14);
ylabel('z [m]','fontsize',14);
legend([h1 h2 h3 h4 h5 h6 h7],'Initial Cw','Cw','Concentration Profile (Observed)','Ct','Cw-ret','Cw-ret+Ct-deg','Ct-deg','Location','southeast');
ylim([-1.0 0]);
xlim([0 0.03]);
set(gca,'fontsize',14,'linewidth',2,'XMinorTick','on','YMinorTick','on');


end

