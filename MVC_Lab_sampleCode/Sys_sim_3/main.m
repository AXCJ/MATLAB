% <<< Linear Quadratic Tracker >>> %
clc; clear ; close all;

Tf = 0.01; % sampling period
Tend = 10; % final simulation time

amp = 1; % Amplitude of square wave
peri = 2; % Period of square wave

% ~~~ system setting
A = [  0  ,   1   ;
     -1.2 , -0.8 ];
 
B = [ 0  ;
      1 ];
  
C = [ 0.2 , 0.5 ];

x_ini = [ 1  ; -1 ];
% ~~~ system setting ~~~ End 

[n,r] = size(B);
p = size(C,1); % the size of the first dimension of C
% n : the number of states  of plant
% r : the number of inputs  of plant
% p : the number of outputs of plant

% ~~~ Controller design
Rc = 1e0*eye(r);
Qc = 1e4*eye(p); % Tracker --> the size of Qc is p by p
[Kc,Pc] = lqr(A, B, C'*Qc*C, Rc); % Tracker --> C'*Q*C
Ec = -inv(Rc)*B'*inv(A - B*Kc)'*C'*Qc;
% ~~~ Controller design ~~~ End

x_c(:,1) = x_ini;
ii = 0; % continuous-time index

% ~~~ Perform simulation
  for t_ = 0 : Tf : Tend
      ii = ii + 1;
      tf(:,ii) = t_; % continuous timespan
    
      % ~~~ Generate a square wave
      cond(:,ii) = tf(:,ii)/peri-floor(tf(:,ii)/peri); % condition
      if (cond(:,ii) < 0.5)
         ref(:, ii) =  amp;
      else
         ref(:, ii) = -amp;    
      end
      % ~~~ Generate a square wave ~~~ End
     
      u_c(:, ii) = Ec*ref(:, ii) - Kc*x_c(:, ii); % Tracker --> the reference input isn't equal to zero
    
      % ~~~ Use R-K's method to solve ODEs ~~~ %
      u(:, ii) = u_c(:, ii);
      options = odeset('RelTol',1e-7);
      tspan = [ t_ , t_+Tf ]; % simulation interval
      [T,X] = ode45(@(tt, xx) plant1(tt, xx, u(:, ii), A, B), tspan, x_c(:,ii), options);
      y_c(:, ii) = C*x_c(:, ii);
    
      if(t_ >= Tend)
          break;
      end
    
      x_c(:,ii+1) = X(end,:)';
    
  end
% ~~~ Perform simulation ~~~ End

% ~~~ Plot fig. ~~~ %
scrsz = get(groot,'ScreenSize');
% [Left bottom width height]
figure('Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])
plot(tf, ref, 'b-', tf, y_c, 'r-.', 'LineWidth', 1.5)
axis([0 Tend -amp*1.1 amp*1.1])
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of \Gamma(t) and y_{c}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf \Gamma','\bf y_{c}')

figure('Position',[scrsz(4)/1.6 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])
plot(tf, x_c(1,:), 'b-', tf, x_c(2,:), 'r-.', 'LineWidth', 1.5)
grid on
title('\bf\fontsize{14}\fontname{Cambria} The time responses of x_{c}(t)')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf x_{c1}','\bf x_{c2}')

