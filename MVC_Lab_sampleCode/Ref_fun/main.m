% <<< Use for-loop to write a reference function >>> %
clc; close all; clear all;
Tf = 0.01; % sampling period
Tend = 10; % final simulation time

ii = 0; % continuous-time index
  for t_ = 0 : Tf : Tend
      ii = ii + 1
      tf(:, ii) = t_; % continuous timespan
      [u1(:, ii) u2(:, ii)] = ref_fun(t_ );
      u(:, ii) = [u1(:, ii) ; u2(:, ii)];
  end

plot(tf, u, 'LineWidth',1.5)
plot(tf, u(1,:), 'b-', tf, u(2,:), 'r-.', 'LineWidth',1.5)
axis([0 Tend -1.2 1.2])
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of control signal')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf u_{1}','\bf u_{2}')
