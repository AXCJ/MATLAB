% Multivariable Control
% 
% HW1
% 
% Do the computer simulation for the following chaotic
% system via two approachs, i.e. (i)simulink and (ii)ode45 respectively.
clc;clear;close all;

%% main switch of figure 1~10
is_fig_1_2_3_4 = 1;
is_fig_5_6_7 = 1;
is_fig_8_9_10 = 1;

%%
axis_fig2_fig4 = [0 60 -30 30 -30 40];

% open('sim_chaotic_system'); % launch the model
% load_system('sim_chaotic_system');
% simMode = get_param('sim_chaotic_system', 'SimulationMode');

Tf = 0.005;
tspan10 = 0:Tf:10;
tspan60 = 0:Tf:60;

x0 = [-10 0 37]; % initial condition

%% sim
%-----config. parameter-----%
% configSet = getActiveConfigSet('sim_chaotic_system');
% configSet.get_param('SolverName')
% get_param(configSet, 'StopTime') % Set StopTime
% set_param('sim_chaotic_system/c','SampleTime', 'Tf') % Set sample time
% set_param('sim_chaotic_system','initialstep','Tf')

% start sim
% set_param('sim_chaotic_system',  'StopTime', '10') % set time stop at 10 second
% [T10_sim,simout_10]=sim('sim_chaotic_system');
[T10_sim,simout_10]=sim('HW1_sim_chaotic_system',tspan10);
% set_param('sim_chaotic_system',  'StopTime', '60') % set time stop at 60 second
% [T60_sim,simout_60]=sim('sim_chaotic_system');
[T60_sim,simout_60]=sim('HW1_sim_chaotic_system',tspan60);
% [tspan_sim,span_simout_10]=sim('sim_chaotic_system',tspan10');

%% ode45
% xdot =chaotic_system(x0);
% [T, dxdt]=ode45(@ (t,y)chaotic_system(tspan10,x0),tspan10,x0);
[tspan10, dxdt_10]=ode45(@ HW1_fun_chaotic_system,tspan10,x0);
[tspan60, dxdt_60]=ode45(@ HW1_fun_chaotic_system,tspan60,x0);


if is_fig_1_2_3_4
%% Figure 1. [Simulink] x1,x2,x3
if 1
figure('Name','[Simulink] x1 , x2 , x3','NumberTitle','on')
subplot(311)
plot(tspan10,simout_10(:,1),'LineWidth',1);
xlabel('t(s)','fontsize',10,'fontweight','b');
title('[Simulink] x1','fontsize',12,'fontweight','b');
% text(t(20,1),x(20,1),'\color{red} \leftarrow \beta',...
%      'HorizontalAlignment','left',...
%      'FontSize',18)
set(gca,'fontsize',10);subplot(312)
subplot(312)
plot(tspan10,simout_10(:,2),'LineWidth',1);
xlabel('t(s)','fontsize',10,'fontweight','b');
title('[Simulink] x2','fontsize',12,'fontweight','b');
subplot(313)
plot(tspan10,simout_10(:,3),'LineWidth',1);
xlabel('t(s)','fontsize',10,'fontweight','b');
title('[Simulink] x3','fontsize',12,'fontweight','b');
end

%% Figure 2. [Simulink] x3-x1-x2
if 1
figure('Name','[Simulink] x3-x1-x2','NumberTitle','on')
plot3(simout_60(:,3),simout_60(:,1),simout_60(:,2));
axis(axis_fig2_fig4)
xlabel('x3','fontsize',10,'fontweight','b');
ylabel('x1','fontsize',10,'fontweight','b');
zlabel('x2','fontsize',10,'fontweight','b');

title('[Simulink] x3-x1-x2','fontsize',12,'fontweight','b');
% close all;
end




%% Figure 3. [ode45] x1 , x2 , x3
if 1
figure('Name','[ode45] x1 , x2 , x3','NumberTitle','on')
subplot(311)
plot(tspan10,dxdt_10(:,1),'LineWidth',1);
xlabel('t(s)','fontsize',10,'fontweight','b');
title('[ode45] x1','fontsize',12,'fontweight','b');
% text(t(20,1),x(20,1),'\color{red} \leftarrow \beta',...
%      'HorizontalAlignment','left',...
%      'FontSize',18)
set(gca,'fontsize',10);subplot(312)
subplot(312)
plot(tspan10,dxdt_10(:,2),'LineWidth',1);
xlabel('t(s)','fontsize',10,'fontweight','b');
title('[ode45] x2','fontsize',12,'fontweight','b');
subplot(313)
plot(tspan10,dxdt_10(:,3),'LineWidth',1);
xlabel('t(s)','fontsize',10,'fontweight','b');
title('[ode45] x3','fontsize',12,'fontweight','b');
end

%% Figure 4. [ode45] x3-x1-x2
if 1
figure('Name','[ode45] x3-x1-x2','NumberTitle','on')
plot3(dxdt_60(:,3),dxdt_60(:,1),dxdt_60(:,2));
axis(axis_fig2_fig4)
xlabel('x3','fontsize',10,'fontweight','b');
ylabel('x1','fontsize',10,'fontweight','b');
zlabel('x2','fontsize',10,'fontweight','b');
title('[ode45] x3-x1-x2','fontsize',12,'fontweight','b');
end
end

% error of sim_x and ode45_x 
err_x= simout_10-dxdt_10;
% abs_err_x = abs(err_x);
% abs_err_x = abs_err_x./(abs(simout_10));
abs_err_x = (abs(err_x))./(abs(simout_10));

if is_fig_5_6_7
%% Figure 5. (|[Simulink] x1 - [ode45] x1|) / |[Simulink] x1|
if 1
figure('Name','(|[Simulink] x1 - [ode45] x1|) / |[Simulink] x1|','NumberTitle','on')
plot(tspan10,abs_err_x(:,1),'LineWidth',1);
xlabel('t(s)','fontsize',10,'fontweight','b');
title('(|[Simulink] x1 - [ode45] x1|) / |[Simulink] x1|','fontsize',12,'fontweight','b');
set(gca,'fontsize',10);
end

%% Figure 6. (|[Simulink] x1 - [ode45] x1|) / |[Simulink] x1|
if 1
figure('Name','(|[Simulink] x2 - [ode45] x2|) / |[Simulink] x2|','NumberTitle','on')
plot(tspan10,abs_err_x(:,2),'LineWidth',1);
xlabel('t(s)','fontsize',10,'fontweight','b');
title('(|[Simulink] x2 - [ode45] x2|) / |[Simulink] x2|','fontsize',12,'fontweight','b');
end

%% Figure 7. (|[Simulink] x1 - [ode45] x1|) / |[Simulink] x1|
if 1
figure('Name','(|[Simulink] x3 - [ode45] x3|) / |[Simulink] x3|','NumberTitle','on')
plot(tspan10,abs_err_x(:,3),'LineWidth',1);
xlabel('t(s)','fontsize',10,'fontweight','b');
title('(|[Simulink] x3 - [ode45] x3|) / |[Simulink] x3|','fontsize',12,'fontweight','b');
end

end

er_x= simout_10-dxdt_10;

if is_fig_8_9_10
%% Figure 8. [Simulink] x1 - [ode45] x1
if 1
figure('Name','[Simulink] x1 - [ode45] x1','NumberTitle','on')
plot(tspan10,simout_10(:,1),'LineWidth',1);
xlabel('t(s)','fontsize',10,'fontweight','b');
title('[Simulink] x1 - [ode45] x1','fontsize',12,'fontweight','b');
set(gca,'fontsize',10);
hold on
plot(tspan10,dxdt_10(:,1),'LineWidth',1);
hold off
end

%% Figure 9. [Simulink] x1 - [ode45] x1
if 1
figure('Name','[Simulink] x1 - [ode45] x1','NumberTitle','on')
plot(tspan10,simout_10(:,2),'LineWidth',1);
xlabel('t(s)','fontsize',10,'fontweight','b');
title('[Simulink] x2 - [ode45] x2','fontsize',12,'fontweight','b');
hold on
plot(tspan10,dxdt_10(:,2),'LineWidth',1);
hold off
end

%% Figure 10. (|[Simulink] x1 - [ode45] x1|) / |[Simulink] x1|
if 1
figure('Name','[Simulink] x1 - [ode45] x1','NumberTitle','on')
plot(tspan10,simout_10(:,3),'LineWidth',1);
xlabel('t(s)','fontsize',10,'fontweight','b');
title('[Simulink] x3 - [ode45] x3','fontsize',12,'fontweight','b');
hold on
plot(tspan10,dxdt_10(:,3),'LineWidth',1);
hold off
end

end