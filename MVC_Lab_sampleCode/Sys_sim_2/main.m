% <<< 各類型輸入與輸出 >>> %
clc; clear all; close all;

Tf = 0.01; % sampling period
Ts = 0.1;  % sampling time
Tend = 10; % final simulation time
sample_div = fix(Ts/Tf);

% ~~~ display
fprintf('sampling period Tf:%2.3f sec\n',Tf);
disp(sprintf('sampling time Ts:%2.3f sec',Ts));
disp(sprintf('sample_div :%2.3f',fix(Ts/Tf)));

% ~~~ system setting
A = [  0  ,   1   ;
     -1.2 , -0.8 ];
 
B = [ 0  ;
      1 ];
  
C = [ 0.2 , 0.5 ];

x_ini = [ 1  ; -1 ];
% ~~~ system setting ~~~ End

x(:, 1) = x_ini;
ii = 0; % continuous-time index
kk = 0; %   discrete-time index

% ~~~ Perform simulation
for t_ = 0 : Tf : Tend
    ii = ii + 1;
    tf(:, ii) = t_; % continuous-time span
    u_c(:, ii) = sin(t_); % 未經取樣之參考輸入 --> 連續訊號
    
    % ~~~ Input sample
    if( mod(ii-1, sample_div) == 0) 
        kk = kk + 1;
        ts(:, kk) = tf(:, ii); %  discrete-time span
        u_d(:, kk) = u_c(:, ii); % 經取樣之參考輸入 --> 離散訊號
    end
    u_zoh(:, ii) = u_d(:, kk); % 先經取樣再 Z.O.H 之參考輸入 --> 連續訊號
    % ~~~ Input sample ~~~ End
    
    u(:, ii) = u_zoh(:, ii);
    options = odeset('RelTol',1e-7, 'AbsTol', 1e-7);
    tspan = [ t_ , t_+Tf ]; % simulation interval
    [T,X] = ode45(@(tt, xx) plant1(tt, xx, u(:, ii), A, B), tspan, x(:, ii), options);
    % Use R-K's method to solve ODEs
    y_c(:, ii) = C*x(:, ii); % 未經取樣之輸出 --> 連續訊號
    
    % ~~~ Output sample
    if( mod(ii-1, sample_div) == 0) 
        y_d(:, kk) = y_c(:, ii); % 經取樣之輸出 --> 離散訊號
    end
    % ~~~ Output sample ~~~ End
    
    y_zoh(:, ii) = y_d(:, kk); % output signal hold
      
    x(:, ii+1) = X(end,:)'; % 下一刻狀態等於 X 之最後一列轉置
end
% ~~~ Perform simulation ~~~ End

% ~~~ Plot fig. ~~~ %
scrsz = get(groot,'ScreenSize');
% [Left bottom width height]
figure('Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])
plot(tf, u_c, '-', ts, u_d, 'o', tf, u_zoh, '-.', 'LineWidth',1.5)
axis([0 Tend -1.1 1.1])
grid on
title('\bf\fontsize{14}\fontname{Cambria} The actuating signal u')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf u_{c}','\bf u_{d}','\bf u_{zoh}')

% figure('Position',[9 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3])
figure('Position',[scrsz(4)/1.6 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])
plot(tf, y_c, '-', ts, y_d, 'o', tf, y_zoh, '-.', 'LineWidth',1.5)
grid on
title('\bf\fontsize{14}\fontname{Cambria} The output y')
xlabel('\bf\fontsize{10} Time (sec.)')
ylabel('\bf\fontsize{10} Amplitude')
legend('\bf y_{c}','\bf y_{d}','\bf y_{zoh}')




