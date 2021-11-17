%% Function to plot soil moisture and concentration conditions in soil matrix and pfd at current time step

function plot_currstate(pfd_n, n_mak,pfd_m, pfd_Ct_deg, pfd_dz, pfd_Ct, pfd_Cw_ret, dz, Ct0,Cm, Ct_deg, Ct, Cw_ret, Cw, dim, istyp, pfd_Cw, pfd_dim, pfd_theta, pfd_z, theta, theta_event, ths, time, z, z_plot)


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

% concentration profile of soil matrix
subplot(3,2,2)
h2 = plot((Cw(1:dim-1).*theta(1:dim-1)),z_plot(1:dim-1)','r-','linewidth',2);
hold on;
h3 = plot((Cw_ret(1:dim-1).*theta(1:dim-1))+(Ct(1:dim-1)./dz),z_plot(1:dim-1)','y-','linewidth',2);
h4 = plot((Cw_ret(1:dim-1).*theta(1:dim-1)),z_plot(1:dim-1)','g-','linewidth',2);
h5 = plot((Cw_ret(1:dim-1).*theta(1:dim-1))+(Ct_deg(1:dim-1)./dz),z_plot(1:dim-1)','c-','linewidth',2,'markersize',2);
h6 = plot(Ct_deg(1:dim-1)./dz,z_plot(1:dim-1)','k--','linewidth',2,'markersize',2);
hold off;
title(['time' num2str(time) ' s' ],'fontsize',14);
xlabel('C [kg/m³]','fontsize',14);
ylabel('z [m]','fontsize',14);
legend([h2 h3 h4 h5 h6],'Cw','Cw-ret+Ct','Cw-ret','Cw-ret+Ct_deg','Ct-deg','Location','southeast');
set(gca,'fontsize',14,'linewidth',2);

%Cm - concentration profile of metabolite in matrix
subplot(3,2,3)
plot(Cm(1:dim-1).*theta(1:dim-1),z_plot(1:dim-1)','y-','linewidth',2);
title(['time' num2str(time) ' s' ],'fontsize',14);
xlabel('Cm','fontsize',14);
ylabel('z [m]','fontsize',14);
set(gca,'fontsize',14,'linewidth',2);

% soil moisture in pfd
subplot(3,2,5)
plot(pfd_theta(1:pfd_dim-1)',pfd_z(1:pfd_dim-1),'r-','linewidth',2);
title(['time' num2str(time) ' s' ],'fontsize',14);
xlabel('Water Content in Pfd','fontsize',14);
ylabel('Pfd_z [m]','fontsize',14);
set(gca,'fontsize',14,'linewidth',2);

% concentration profile of pfd
subplot(3,2,6)
h2 = plot((((((pfd_n(1) * n_mak) .* pfd_theta(1:pfd_dim-1)) .* pfd_m) ./ 1000) .* pfd_Cw(1:pfd_dim-1)) ./ pfd_dz,pfd_z(1:pfd_dim-1),'r-','linewidth',2);
hold on;
h3 = plot(((((((pfd_n(1) * n_mak) .* pfd_theta(1:pfd_dim-1)) .* pfd_m) ./ 1000) .* pfd_Cw_ret(1:pfd_dim-1)) ./ pfd_dz) + (pfd_Ct(1:pfd_dim-1) ./ pfd_dz),pfd_z(1:pfd_dim-1)','y-','linewidth',2);
h4 = plot((((((pfd_n(1) * n_mak) .* pfd_theta(1:pfd_dim-1)) .* pfd_m) ./ 1000) .* pfd_Cw_ret(1:pfd_dim-1)) ./ pfd_dz, pfd_z(1:pfd_dim-1)','g-','linewidth',2);
h5 = plot(((((((pfd_n(1) * n_mak) .* pfd_theta(1:pfd_dim-1)) .* pfd_m) ./ 1000) .* pfd_Cw_ret(1:pfd_dim-1)) ./ pfd_dz) + (pfd_Ct_deg(1:pfd_dim-1) ./ pfd_dz),pfd_z(1:pfd_dim-1)','c-','linewidth',2,'markersize',2);
h6 = plot(pfd_Ct_deg(1:pfd_dim-1) ./ pfd_dz,pfd_z(1:pfd_dim-1)','k--','linewidth',2,'markersize',2);
hold off;
title(['time' num2str(time) ' s' ],'fontsize',14);
xlabel('soil concentration [kg/m³]','fontsize',14);
ylabel('z [m]','fontsize',14);
% legend([h2 h3 h4 h5 h6],'pfdCw','pfdCw-ret+pfdCt','pfdCw-ret','pfdCw-ret+pfdCtdeg','pfdCt-deg','Location','southeast');
% xlim([0 2 ]);
set(gca,'fontsize',14,'linewidth',2);

end

