%% Function to plot the final concentration profile

function plot_age(dim, mtx_avg_age, z_plot)

figure;

hold on;
h1=plot(mtx_avg_age(1:dim-1),z_plot(1:dim-1)','k-','linewidth',2,'markersize',2);
hold off;
title('Age Profile at End of Simulation','fontsize',14);
xlabel('Age [s]','fontsize',14);
ylabel('z [m]','fontsize',14);
legend(h1,'Age','Location','southeast');
ylim([-1.0 0]);
% xlim([0 0.03]);
set(gca,'fontsize',14,'linewidth',2,'XMinorTick','on','YMinorTick','on');


end

