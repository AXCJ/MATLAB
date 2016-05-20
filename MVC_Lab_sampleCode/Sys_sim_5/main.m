% <<< Continuous-time system plus observer simulation >>> %
clc; clear all; close all;

Tf = 0.01; % sampling period
Tend = 10; % final simulation time

% ~~~ System setting
A = [  0  ,   1   ;
     -1.2 , -0.8 ];
 
B = [ 0  ;
      1 ];
  
C = [ 0.2 , 0.5 ];

x_ini = [ 1  ; -1 ];
% ~~~ System setting ~~~ End  

[n,r] = size(B);
p = size(C,1); % the size of the first dimension of C 
% n : the number of states  of plant
% r : the number of inputs  of plant
% p : the number of outputs of plant

% ~~~ Observer setting
Ao = [  0  ,   1   ;
      -1.2 , -0.8 ];
 
Bo = [ 0  ;
       1 ];
  
Co = [ 0.2 , 0.5 ];

xo_ini = [ 0  ; 0 ];
% ~~~ Observer setting ~~~ End  

[no,ro] = size(Bo);
po = size(Co,1);
% no : the number of states  of observer
% ro : the number of inputs  of observer
% po : the number of outputs of observer

% ~~~ Observer design
Ro = 1e0*eye(ro);
Qo = 1e4*eye(no); % Regulator --> the size of Qo is no by no

[Ko,Po] = lqr(A', C', Qo, Ro); % Regulator --> Qo
Ko = Ko'; % (A, B ,Kc) <---> (A', C', Ko') for continuous-time system
% Design Ko such that ( A - Ko*C ) is asymptotically stable
% ~~~ Observer design ~~~ End

x_c(:,1) = x_ini;
x_o(:,1) = xo_ini;
ii = 0;

% ~~~ Perform simulation
  for t_ = 0 : Tf : Tend
      ii = ii + 1
      tf(:,ii) = t_; % continuous timespan
      u(:, ii) = ref_fun(t_);
      options = odeset('RelTol',1e-7);
      tspan = [ t_ , t_+Tf ];
      
      [Tc,Xc] = ode45(@(tt, xx) plant1(tt, xx, u(:, ii), A, B), tspan, x_c(:,ii), options);
      y_c(:, ii) = C*x_c(:, ii);
    
      [To,Xo] = ode45(@(tt, xx) observer1(tt, xx, u(:, ii), Ao, Bo, Co, Ko, y_c(:, ii)), tspan, x_o(:,ii), options);
      y_o(:, ii) = Co*x_o(:, ii);
     
      if(t_ >= Tend)
          break;
      end
    
      x_c(:,ii+1) = Xc(end,:)';
      x_o(:,ii+1) = Xo(end,:)';
    
  end
% ~~~ Perform simulation ~~~ End

% ~~~ Plot fig. ~~~ %
figure(1)
plot(tf, u, '-', tf, y_c, ':', tf, y_o,'-.', 'LineWidth',1.5)
axis([0 Tend -0.8 2.2])
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of u(t) , y_{c}(t) and y_{o}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf u','\bf y_{c}','\bf y_{o}')

figure(2)
plot(tf, x_c(1,:), 'b-', tf, x_o(1,:),'r-.', 'LineWidth',1.5)
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of x_{1}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf x_{c1}','\bf x_{o1}')

figure(3)
plot(tf, x_c(2,:), 'b-', tf, x_o(2,:),'r-.', 'LineWidth',1.5)
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of x_{2}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf x_{c2}','\bf x_{o2}')
