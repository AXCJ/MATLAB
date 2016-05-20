clear all; clc; close all;
Tf = 1e-1;
Tend = 100;
sigma0 = 1e5;
sigma1 = sqrt(1e5);
sigma2 = 0.4;
Fc = 1;
Fs = 1.5;
vs = 0.001;
m = 1;
Kv = 6;
Kp = 3;
Ki = 4;

options = simset('solver','ode23s');
[t,~,pos] = sim('model',(0:Tf:Tend), options);
figure
h=plot(t,pos, t,ones(1,length(t)));
axis([-inf inf -inf inf])
set(gca, 'Box', 'off')
xlabel('\bf\fontsize{20} Time (sec.)')
ylabel('\bf\fontsize{20} position(m)')
% set_param('simulink','Lock','off')

