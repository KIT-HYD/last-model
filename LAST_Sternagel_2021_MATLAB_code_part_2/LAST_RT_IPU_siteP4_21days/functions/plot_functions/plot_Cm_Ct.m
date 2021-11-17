%% Function to plot the concentration profile of the metabolite and solute in adsorbing phase

function plot_Cm_Ct(Ct0, Cm, theta, z_plot, dim, dz)

figure;

hold on;
h2=plot(Cm(1:dim-1) ./ dz,z_plot(1:dim-1)','k-','linewidth',2,'markersize',2);

h5 = plot(Ct0(1:dim-1).*theta(1:dim-1),z_plot(1:dim-1)','b-','linewidth',2,'markersize',2);
hold off;
title('Solute Concentration','fontsize',14);
xlabel('Concentration [kg/m³]','fontsize',14);
ylabel('z [m]','fontsize',14);
legend([h2 h5],'Cm','Ct0','Location','southeast');
ylim([-1.0 0]);
% xlim([0 0.03]);
set(gca,'fontsize',14,'linewidth',2,'XMinorTick','on','YMinorTick','on');

end

