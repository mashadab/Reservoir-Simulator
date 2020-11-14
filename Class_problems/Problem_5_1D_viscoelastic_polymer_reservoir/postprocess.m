function postprocess(P_plot,visc_plot,vel_plot, shear_plot,numerical, reservoir, fluid)

tplot = [1;2;3;4;5;7;10;12;15; 20; 25; 35; 50]; %; 75; 87; 100; numerical.t_final-5];   %days

tplot = [2;3; 5;10];
%colors = {'k-o'; 'k-o'; 'k-o'; 'k-o';'k-o'; 'k-o'; 'k-o'; 'k-o'; 'k-o'; 'k-o'; 'k-o'; 'k-o'};

% for i=1:length(tplot)
%     hold on
%     m = round(tplot(i)/numerical.dt);
%     plot(numerical.x, P_plot(:,m),'k-')
% end

figure(2)
plot(numerical.x, P_plot(:,2),'k-',numerical.x, P_plot(:,10),'k-',numerical.x, P_plot(:,50),'k-', numerical.x, P_plot(:,100),'k-')
set(gca,'fontsize',13)
xlabel ('x (feet)')
ylabel ('Pressure (psi)')
legend ('2 days', '10 days','30 days','0')
whitebg('w')
print -djpeg -r300 'Pressure.jpg'
figure(3)

plot(numerical.x, vel_plot(:,2),'k-',numerical.x, vel_plot(:,10),'k-',numerical.x, vel_plot(:,50),'k-', numerical.x, vel_plot(:,100),'k-')
set(gca,'fontsize',13)
xlabel ('x (feet)')
ylabel ('Velocity (ft/day)')
%ylim([0 100])
legend ('2 days', '10 days','30 days','0')
whitebg('w')
print -djpeg -r300 'Velcoity.jpg'
figure(4)

plot(numerical.x, shear_plot(:,2),'k-',numerical.x, shear_plot(:,10),'k-',numerical.x, shear_plot(:,50),'k-', numerical.x, shear_plot(:,100),'k-')
set(gca,'fontsize',13)
xlabel ('x (feet)')
ylabel ('Shear Rate (1/s)')
%ylim([0 100])
legend ('2 days', '10 days','30 days','0')
whitebg('w')
print -djpeg -r300 'Shear_rate.jpg'
figure(5)

% for i=1:length(tplot)
%     hold on
%     m = round(tplot(i)/numerical.dt);
%     plot(numerical.x, visc_plot(:,m),'k-')
%     %plot(numerical.x(1:numerical.N-1), velocity_half,'k-')
%     %plot(numerical.x(1:numerical.N-1), shear_half,'k-')
% end
plot(numerical.x, visc_plot(:,2),'k-',numerical.x, visc_plot(:,10),'k-',numerical.x, visc_plot(:,50),'k-', numerical.x, visc_plot(:,100),'k-')
set(gca,'fontsize',13)
xlabel ('x (feet)')
ylabel ('Viscosity (cp)')
%ylim([0 100])
legend ('2 days', '10 days','30 days','0')
whitebg('w')
print -djpeg -r300 'Viscosity.jpg'
%figure(3)

% %plot(x, P_plot(:,n1),'-o', x,P_plot(:,n2),'-o', x, P_plot(:,n3), '-o',x,P_plot(:,n4),'-o')
% set(gca,'fontsize',13)
% xlabel ('x (feet)')
% ylabel ('Pressure (psi)')
% %ylabel ('Apparent viscosity (cp)')
% %ylabel ('Velocity (ft/day)')
% %ylabel ('Shear rate (1/s)')
% legend ('2 days', '10 days','30 days','0')
% whitebg('w')
% print -djpeg -r300 'Viscosity(Jacobian).jpg'