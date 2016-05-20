
clc;clear all;close all;
Linewidth = 1;
scrsz = get(groot,'ScreenSize');
Tf = 5e-3; % sampling period
Tend = 3; % final simulation time
Ts = 0.01; % sampling time
sample_div = fix(Ts/Tf);

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
x_ini = [ 0.3 -0.2 0 -0.1 0.2]';
% ~~~ system setting ~~~ End


[n,r] = size(B);
p = size(C,1); % the size of the first dimension of C
% n : the number of states of plant
% r : the number of inputs of plant
% p : the number of outputs of plant

% ~~~ Controller design
Rc = 1e0*eye(r);
Qc = 1e6*eye(p); % Tracker --> the size of Qc is p by p
[Kc,Pc] = lqr(A, B, C'*Qc*C, Rc); % Tracker --> C'*Q*C
Ec = -inv(Rc)*B'*inv(A - B*Kc)'*C'*Qc;
% ~~~ Controller design ~~~ End

xc(:,1) = x_ini;
ii = 0; % continuous-time index
kk = 0; %   discrete-time index
% ~~~ Perform simulation
for t_ = 0 : Tf : Tend
    ii = ii + 1;
    ref(: , ii) = ref_fun(t_)'; % get reference input
    ref_star(: , ii) = ref_fun(t_+Tf)'; % ref_star(:, ii) is equal to ref(: , ii+1)

%    fprintf('t_= %.7f\n',t_)
    tf(: , ii) = t_; % continuous-time span


    uc(: , ii) = Ec*ref_star(: , ii) - Kc*xc(: , ii); % 未經取樣之輸出 --> 連續訊號
    
    options = odeset('RelTol',1e-7);
    tspan = [ t_ , t_+Tf ]; % simulation interval
    % Use R-K's method to solve ODEs
    [T,X] = ode45(@(tt, xx) plant1(tt, xx, uc(:, ii), A, B), tspan, xc(:, ii), options);
    yc(:, ii) = C*xc(:, ii); % 未經取樣之輸出 --> 連續訊號
    if ( mod(ii-1, sample_div) == 0 )
        kk = kk + 1;
        ts(:, kk) = tf(:, ii); %  discrete timespan
        xc_s(:, kk) = xc(:, ii);
        ref_s(:, kk) = ref(:, ii);
        uc_s(:, kk) = uc(:, ii);
        yc_s(:, kk) = yc(:, ii);
    end
    if(t_ >=Tend)
        break;
    end
    
    xc(:, ii+1) = X(end,:)'; % 下一刻狀態等於 X 之最後一列轉置
end


% ~~~ Plot fig. ~~~ %
% [Left bottom width height]

figure('Name', 'ref1 vs. yc1' ,'Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Left_Top of screen
plot(ts , ref_s(1 , :) , ts , yc_s(1 , :) , 'LineWidth' , Linewidth)
legend({'\bfu_{c1}','\bfy_{c1}'} , 'fontsize',15)
xlabel('t(s)','fontsize',20)
title('y_{c1} vs. r_{1}','fontsize',30)
hh=gca; set(hh,'FontSize',20)

figure('Name', 'ref2 vs. yc2' ,'Position',[scrsz(3)/3 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Middle_Top of screen
plot(ts , ref_s(2 , :), ts , yc_s(2 , :) , 'LineWidth' , Linewidth)
legend({'\bfu_{c2}','\bfy_{c2}'} , 'fontsize',15)
xlabel('t(s)','fontsize',20)
title('y_{c2} vs. r_{2}','fontsize',30)
hh=gca; set(hh,'FontSize',20)

figure('Name', 'u1' ,'Position',[9 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3]) % Left_Bottom of screen
h111=plot(ts , uc_s(1 , :))
hh=gca; set(hh,'FontSize',20)
xlabel('t(s)','fontsize',20)
title('u_{1}','fontsize',30)
% set(h111,'LineWidth',10)

figure('Name', 'u2' ,'Position',[scrsz(3)/3 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3]) % Middle_Bottom of screen
plot(ts , uc_s(2 , :))
hh=gca; set(hh,'FontSize',20)
xlabel('t(s)','fontsize',20)
title('u_{2}','fontsize',30)

sys = ss(A,B,C,0); % open-loop system
figure('Name', 'Open-loop Pole-Zero Map' ,'Position',[scrsz(3)/1.5 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Right_Top of screen
h1 = pzplot(sys);
h1o = getoptions(h1); % Get options for plot.
h1o.Title.FontSize = 18;
h1o.Title.String = 'Open-loop Pole-Zero Map';
h1o.XLabel.FontSize = 15;
h1o.YLabel.FontSize = 20;
h1o.TickLabel.FontSize = 15;
setoptions(h1,h1o); % Apply options to plot.

syscl = ss(A-B*Kc,B*Ec,C,0); % closed-loop system
figure('Name', 'Closed-loop Pole-Zero Map' ,'Position',[scrsz(3)/1.5 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3]) % Right_Bottom of screen
h2 = pzplot(syscl);
h2o = getoptions(h2); % Get options for plot.
h2o.Title.FontSize = 20;
h2o.Title.String = 'Closed-loop Pole-Zero Map';
h2o.XLabel.FontSize = 15;
h2o.YLabel.FontSize = 20;
h2o.TickLabel.FontSize = 15;
setoptions(h2,h2o); % Apply options to plot.



%%%~~~ another solution ~~~%%%

% Tf = 1e-2;
% Td_seq = 0 : Tf : Tend;
% 
% figure('Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Left_Top of screen
% options = odeset('RelTol',1e-7);
% [T,x_dott]=ode45(@(t,x) plant2(t,x,Kc,Ec,A,B,ref) , Td_seq , x_ini , options);
% Y_ode=C *x_dott';


% figure('Position',[scrsz(3)/3 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Middle_Top of screen
% plot(Td_seq , ref(1 , 1:301) , Td_seq , Y_ode(1 , :) , '-.' , 'LineWidth' ,1 )
% legend('r_{1}','y_{c1}')
% xlabel('t(s)')
% title('y_{c1} vs. r_{1}')
% % ylabel('y_{c1} vs. r_{c1}')
% % set(get(gca,'YLabel'),'Rotation',0);
% % set(get(gca,'YLabel'),'Position',get(get(gca,'YLabel'),'Position') - [0.090 0.01 0])
% 
% figure('Position',[scrsz(3)/1.5 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Right_Top of screen
% plot(Td_seq , ref(2 , 1:301), Td_seq , Y_ode(2 , :) , '-.' , 'LineWidth' , 1)
% legend('r_{2}','y_{c2}')
% xlabel('t(s)')
% title('y_{c2} vs. r_{2}')
% % ylabel('y_{c2} vs. r_{c2}')
% % set(get(gca,'YLabel'),'Rotation',0);
% % set(get(gca,'YLabel'),'Position',get(get(gca,'YLabel'),'Position') - [0.090 -.01 0])





