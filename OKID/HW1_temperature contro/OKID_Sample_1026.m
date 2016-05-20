clc; clear; close all;

load('OKIDData.mat')
% load SystemData.mat % load original system
sampleData = 1000; % number of sample data 
deletePoint = 10; % from the first data to the data you choose won't be shown in figure
showPoint = 100; % total plot point
plotPoint = deletePoint + showPoint;

MarkovOrder = 8;

x_initial = [-351.2222;11.2992;8.3836;8.0217;-4.4119;7.6555;37.6572;-23.7958];

u_star = [5000;3.4321;1090;24;2.1837;130;4.4834;380];

Kss = [-0.0722;0.0896;-0.0220;-0.0220;-0.0234;-0.0257;-0.0360;0.0164];

input_weight = [1 30 1 5 40 5 15 2]; % input wighting make input to be same level

kk = 0;
% input generator
for t = 1:1:sampleData
    kk = kk + 1;
    u_randn = [    input_weight(1)*randn(1,1);
                         input_weight(2)*randn(1,1);
                         input_weight(3)*randn(1,1);
                         input_weight(4)*randn(1,1);
                         input_weight(5)*randn(1,1);
                         input_weight(6)*randn(1,1);
                         input_weight(7)*randn(1,1);
                         input_weight(8)*randn(1,1)];
    esys(:,kk) = -0.0073 + 13.3067*randn(1,1);

    u(:,kk) = u_star + u_randn;
end

% original output generator
kk = 0;
xd(:,1) = x_initial;
for t = 1:1:sampleData
    kk = kk + 1;
    xd(:,kk+1) = G*xd(:,kk) + H*u(:,kk) + Kss*esys(:,kk);
    yd(:,kk) = C*xd(:,kk);
end


Ud = u;
Yd = yd;

% Ts = 0.1;
% Tend = 10;
% Sam_t = 0:Ts:Tend-Ts;
% T_length = length(Sam_t);

% OKID the system
IDset_ =  struct('MarkovOrder', [MarkovOrder], 'Alpha', [20], 'Beta', [20], 'n', 0.002, 'MinRA', 'era');
[G_ok,H_ok,C_ok,D_ok,F_ok,Sigma, er, M, Ob_CanF] = ...
    OKID_fun_WT02_6(Ud, Yd, IDset_, []);

% OKID_system output generator 
x_o(:,1) = pinv(C_ok)*(C*x_initial);
for kk = 1 : sampleData
    y_o(:,kk) = C_ok*x_o(:,kk);
    x_o(:,kk+1) = G_ok*x_o(:,kk)  + H_ok*Ud(:,kk) + F_ok*(y_o(:,kk) - Yd(:,kk)) ;
end
% error between original & OKIDed system.
error_Y = Yd - y_o;

ii = 0;
for t = 0 : 1 : sampleData-1
    ii = ii + 1;
    td(:,ii) = ii;
end

% save OKID_esys_without.mat G_ok H_ok C_ok D_ok F_ok

scrsz = get(groot,'ScreenSize');
% [Left bottom width height]
figure('Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Left_Top of screen
plot(td(deletePoint:plotPoint),error_Y(deletePoint:plotPoint)) 
figure('Position',[scrsz(3)/3 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Middle_Top of screen
plot(td(deletePoint:plotPoint),Yd(deletePoint:plotPoint),td(deletePoint:plotPoint),y_o(deletePoint:plotPoint),'.-')
legend('y_{d}','y_{o}')

% figure('Position',[scrsz(3)/1.5 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3]) % Right_Bottom of screen
figure('Position',[scrsz(3)/1.5 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Right_Bottom of screen
plot(td(deletePoint:plotPoint),error_Y(deletePoint:plotPoint)) 
figure('Position',[9 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3]) % Left_Bottom of screen
plot(td(deletePoint:plotPoint),Yd(deletePoint:plotPoint),td(deletePoint:plotPoint),y_o(deletePoint:plotPoint),'.-') 
legend('y_{d}','y_{o}')
eig(G)
eig(G_ok)
% pause
% clc; clear; close all;
% end