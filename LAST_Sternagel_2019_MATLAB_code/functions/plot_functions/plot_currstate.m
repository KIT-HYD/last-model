%% Function to plot soil moisture and concentration conditions in soil matrix and pfd at current time step

function plot_currstate(Cw, dim, istyp, pfd_Cw, pfd_dim, pfd_theta, pfd_z, theta, theta_event, ths, time, z)


% event and pre-event water profile of soil matrix
subplot(3,2,1);
% h1=plot(theta_event(1:dim-1)',z(1:dim-1),'g-','linewidth',2);
h2=plot(theta(1:dim-1)',z(1:dim-1),'r-','linewidth',2);
title(['time' num2str(time) ' s' ],'fontsize',14);
xlabel('theta [-]','fontsize',14);
ylabel('z [m]','fontsize',14);
% legend([h1 h2],'Event Water)', 'Pre-Event Water','Location','southeast');
xlim([0 ths(istyp(1)) ]);
set(gca,'fontsize',14,'linewidth',2);

% total water profile of soil matrix
% subplot(3,2,2)
% h2=plot(theta(1:dim-1)'+theta_event(1:dim-1)',z(1:dim-1),'r-','linewidth',2);
% title(['time' num2str(time) ' s' ],'fontsize',14);
% xlabel('theta [-]','fontsize',14);
% ylabel('z [m]','fontsize',14);
% legend(h2,'Pre-/Event Water(total)','Location','southeast');
% xlim([0 ths(istyp(1)) ]);
% set(gca,'fontsize',14,'linewidth',2);

% concentration profile of soil matrix
subplot(3,2,3)
h2=plot(Cw(1:dim-1)',z(1:dim-1),'r-','linewidth',2);
title(['time' num2str(time) ' s' ],'fontsize',14);
xlabel('Concentration in Soil Matrix [kg/m³]','fontsize',14);
ylabel('z [m]','fontsize',14);
legend(h2,'Solute Concentration','Location','southeast');
set(gca,'fontsize',14,'linewidth',2);

% soil moisture in pfd
subplot(3,2,4)
plot(pfd_theta(1:pfd_dim-1)',pfd_z(1:pfd_dim-1),'r-','linewidth',2);
title(['time' num2str(time) ' s' ],'fontsize',14);
xlabel('Water Content in Pfd','fontsize',14);
ylabel('Pfd_z [m]','fontsize',14);
set(gca,'fontsize',14,'linewidth',2);

% concentration profile of pfd
subplot(3,2,5)
h2=plot(pfd_Cw(1:pfd_dim-1)',pfd_z(1:pfd_dim-1),'r-','linewidth',2);
title(['time' num2str(time) ' s' ],'fontsize',14);
xlabel('Concentration in Pfd [kg/m³]','fontsize',14);
ylabel('z [m]','fontsize',14);
legend(h2,'Solute Concentration','Location','southeast');
xlim([0 2 ]);
set(gca,'fontsize',14,'linewidth',2);

end

