
clc;clear;close all;

scrsz = get(groot,'ScreenSize');

Tf = 1e-2; % sampling period
Tend = 3; % final simulation time
Ts = 2e-1; % sampling time
sample_div = fix(Ts/Tf);
Td_seq = 0 : Tf : Tend;

% ~~~ display
fprintf('sampling period Tf:%2.3f sec\n',Tf);
disp(sprintf('sampling time Ts:%2.3f sec',Ts));
disp(sprintf('sample_div :%2.3f',fix(Ts/Tf)));

% ~~~ system setting
A = [ 0.9   -2.1     0.3     0.5     0.9;
      0.7    0.2     0.3       0     0.7;
     -1.3    0.5    -1.1    -2.3    -0.2;
     -0.3    0.8     1.7    -1.2    -0.3;
     -3.5   -4.3    -0.7       0    -8.3];
B = [   1   -0.4;
     -1.7   -1.7;
     -0.2    1.2;
      0.7    0.1;
      0.9    1.4];
C = [1  -2  -1   0   3;
     0   2   1   0   4];
x_ini = [ 0.3;
             -0.2;
                0;
             -0.1;
              0.2];
% ~~~ system setting ~~~ End


% n : the number of states  of plant
% r : the number of inputs  of plant
[n,r] = size(B);
% p : the number of outputs of plant
p = size(C,1); % the size of the first dimension of C


% ~~~ Continous-time tracker design
R = eye(r);
Q = 5e6*eye(p); % Tracker --> the size of Q is p by p
[Kc,P,e] = lqr(A,B,C'*Q*C,R,0); % Tracker --> C'*Qc*C
Ec = -inv(R)*B'*inv((A - B*Kc)')*C'*Q;% Ec = -R\B'/((A-B*Kc)')*C'*Q;
% ~~~ Continous-time tracker design ~~~ End
% ~~~ Discrete-time tracker design
[G,H] = c2d(A,B,Ts);
Kd = inv(eye(2)+Kc*H)*Kc*G;
Ed = inv(eye(r) + Kc*H)*Ec;
% ~~~ Discrete-time tracker design ~~~ End



% ~~~ Perform Continuous-time simulation
ii = 0; % continuous-time index
% ~~~ get reference input
for t_ = 0 : Tf : Tend+Tf
    ii = ii + 1;
    ref(: , ii) = ref_fun(t_)'; 
end
options = odeset('RelTol',1e-7); 
[Tc,Xc]=ode45(@(t,x) plant2(t,x,Kc,Ec,A,B,ref) , Td_seq , x_ini , options);
Yc = C *Xc';
% ~~~perform continuous time system ~~~End
figure('Position',[scrsz(3)/3 scrsz(4)/6 scrsz(3)/3 scrsz(4)/3])
plot(Td_seq , ref(1 , 1:301) , Td_seq , Yc(1 , :) , '-.' , 'LineWidth' ,1 )
legend('r_{1}','y_{c1}')
xlabel('t(s)')
title('y_{c1} vs. r_{1}')
Xc = Xc';
figure('Position',[8 scrsz(4)/6 scrsz(3)/3 scrsz(4)/3])
plot(Td_seq , Xc(:,:) , '-.' , 'LineWidth' ,1 )
legend('r_{1}','y_{c1}')
xlabel('t(s)')
title('y_{c1} vs. r_{1}')



ii = 0; % continuous-time index
kk = 0; % sampling time index
% ~~~ Perform Discrete-time simulation
for t_ = 0 : Tf : Tend
    ii = ii + 1;

    tf(:,ii) = t_; % continuous-time span


%     uc(:,ii) = Ec*ref(:,ii) - Kc*xc(:,ii);
    % ~~~ u sample
    if( mod(ii-1, sample_div) == 0) 
        kk = kk+1;
        
        ts(:,kk) = tf(:,ii); % sampling time span
        uk(:,kk) = Ec*ref(:,ii+1) - Kc*Xc(:,ii); % 經取樣之參考輸入 --> 離散訊號 (下一點r)
        xd(:,ii) = Xc(:,ii);
    end
    uk_zoh(:,ii) = uk(:,kk); % 先經取樣再 Z.O.H 之參考輸入 --> 連續訊號
    % ~~~ u sample ~~~ End
    
    options = odeset('RelTol',1e-7, 'AbsTol', 1e-7);
    tspan = [ t_ , t_+Tf ]; % simulation interval
    [T,X] = ode45(@(tt, xx) plant1(tt, xx, uk_zoh(:, ii), A, B), tspan, xd(:, ii), options); % X_ini 持續更新
    % Use R-K's method to solve ODEs
    yd(:, ii) = C*xd(:, ii); % 未經取樣之輸出 --> 連續訊號
    
    % ~~~ Output sample
%     if( mod(ii-1, sample_div) == 0) 
%         yk(:,kk) = yc(:,ii);  % 經取樣之輸出 --> 離散訊號
%     end
    % ~~~ Output sample ~~~ End
%     y_zoh(:,ii) = yk(:,kk); % output signal hold
    if t_ >= Tend
        break;
    end
    
    xd(: , ii+1) = X(end)'; % 下一刻狀態等於 X 之最後一列轉置
end

ref(: , end) = [];
% ~~~ Plot fig. ~~~ %

scrsz = get(groot,'ScreenSize');
% [Left bottom width height]
figure('Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])
plot(tf,yd(1,:),tf,Yc(1,:))
legend('y_{d1}','Y_{c1}')
figure('Position',[scrsz(3)/3 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])
plot(tf,yd(1,:),tf,Yc(1,:))
legend('y_{d2}','Y_{c2}')


% figure('Position',[9 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3])
% figure('Position',[scrsz(3)/1.5 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3])
% plot(tf, yc, '-', ts, yk, '.', tf, y_zoh, '-.')
% title('\bf\fontsize{14}\fontname{Cambria} The output y')
% xlabel('\bf\fontsize{10} Time (sec.)')
% ylabel('\bf\fontsize{10} Amplitude')
% legend('\bf y_{c}','\bf y_{k}','\bf y_{zoh}')
















