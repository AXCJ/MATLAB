   % <<< Use R-K method and for-loop to solve ODEs >>> %
clc; clear all; close all;

Tf = 0.01; % sampling period
Tend = 20; % final simulation time

% ~~~ system setting
A = [  0  ,   1   ;
     -1.2 , -0.8 ];
 
B = [ 0  ;
      1 ];
  
C = [ 0.2 , 0.5 ];

x_ini = [ 1  ;
         -1 ];
% ~~~ system setting ~~~ End 

% ~~~ Perform simulation
x(:, 1) = x_ini; 
ii = 0; % continuous-time index
  for t_ = 0 : Tf : Tend
      ii = ii + 1;
      tf(:, ii) = t_; % timespan
      u(:, ii) = 1; % step input
      options = odeset('RelTol',1e-7, 'AbsTol', 1e-7);
      tspan = [ t_ , t_+Tf ]; % simulation interval
      [T, X] = ode45(@(tt, xx) plant1(tt, xx, u(:, ii), A, B), tspan, x(:, ii), options);
      y(:, ii) = C*x(:, ii);
      
      if(t_ >= Tend)
          break;
      end
      
      x(:, ii+1) = X(end,:)'; % next state
    
  end
tspan = 0 : Tf : Tend;
xdot = @(t,x) A*x + B*u(floor(tspan/tf)+1);
[T, X] = ode45(xdot, tspan , x_ini , options);
error_x = x-X';
plot(tf,error_x)
figure
% ~~~ Perform simulation ~~~ End

% --- Plot fig. --- %
plot(tf, y, 'LineWidth',1.5)
grid on
title('\bf\fontsize{14}\fontname{Cambria} The step response')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
