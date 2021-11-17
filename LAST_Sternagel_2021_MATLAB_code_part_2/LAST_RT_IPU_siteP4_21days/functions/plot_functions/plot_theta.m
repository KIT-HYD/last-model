%% Function to plot the final soil moisture profile and additional parameters (v, Pec, D)

function plot_theta(D, dim, dz, theta, theta_init, theta_event, time, v, z)

figure;

% soil mositure profile, separated in the proportions of event and pre-event water
subplot(2,2,1);
hold on;
h1 = plot(theta(1:dim-1)'+theta_event(1:dim-1)',z(1:dim-1),'r-','linewidth',2); % well-mixed water plus the proportion of event water what is still not mixed
h2 = plot(theta(1:dim-1)',z(1:dim-1),'y-','linewidth',2); % only well-mixed water
h3 = plot(theta_init(1:dim-1),z(1:dim-1),'k+','linewidth',2,'markersize',2); % initial soil moisture profile
hold off;
title(['time' num2str(time) ' s' ],'fontsize',14);
xlabel('theta_init [-]','fontsize',14);
ylabel('z [m]','fontsize',14);
legend([h1 h2 h3],'Mixed+Event Water','Well-Mixed Water','Initial Moisture','Location','southeast');
ylim([-1.0 0 ]);
set(gca,'fontsize',14,'linewidth',2,'XMinorTick','on','YMinorTick','on');

% drift velocity
subplot(2,2,2)
h2=plot(smooth(v,4),z(1:dim),'r-','linewidth',2);
xlabel('velocity [m/s] ','fontsize',14);
ylabel('z [m]','fontsize',14);
legend(h2,'Drift Velocity','Location','northeast');
set(gca,'fontsize',14,'linewidth',2)
ylim([z(length(z)) 0 ]);

% grid peclet number
subplot(2,2,3)
plot(smooth(v*min(abs(dz))./D,4),z(1:dim),'b-','linewidth',2);
title(['Peclet number ' num2str(max(v)*max(abs(z))./max(D))],'fontsize',14);
xlabel('Grid Peclet Number ','fontsize',14);
ylabel('z [m]','fontsize',14);
set(gca,'fontsize',14,'linewidth',2)
ylim([z(length(z)) 0 ]);

% diffussion coefficient
subplot(2,2,4)
plot(smooth(D,4),z(1:dim),'b-','linewidth',2);
xlabel('Diffusion Coeff. [m^2/s] ','fontsize',14);
ylabel('z [m]','fontsize',14);
set(gca,'fontsize',14,'linewidth',2);
ylim([z(length(z)) 0 ]);


end

