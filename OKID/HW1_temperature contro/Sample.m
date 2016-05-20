clc; clear; close all;

load OKIDData

x_initial = [-351.2222;11.2992;8.3836;8.0217;-4.4119;7.6555;37.6572;-23.7958];

u_star = [5000;3.4321;1090;24;2.1837;130;4.4834;380];
u_star_3percent = u_star*0.03;
Kss = [-0.0722;0.0896;-0.0220;-0.0220;-0.0234;-0.0257;-0.0360;0.0164];

kk = 0;
for t = 1:1:100
    kk = kk + 1;
    u_randn = [   u_star_3percent(1)*randn(1,1);
                  u_star_3percent(2)*randn(1,1); 
                  u_star_3percent(3)*randn(1,1);
                  u_star_3percent(4)*randn(1,1);
                  u_star_3percent(5)*randn(1,1);
                  u_star_3percent(6)*randn(1,1);
                  u_star_3percent(7)*randn(1,1);
                  u_star_3percent(8)*randn(1,1)];
    esys(:,kk) = -0.0073 + 13.3067*randn(1,1);
%     13.3067

    u(:,kk) = u_star + u_randn; % input
end

kk = 0;
xd(:,1) = x_initial;
for t = 1:1:100
    kk = kk + 1;
    % simulink 
    % x with disturbance
    xd(:,kk+1) = G*xd(:,kk) + H*u(:,kk) + Kss*esys(:,kk);
    % y without disturbance
    yd(:,kk) = C*xd(:,kk); 
end
Ud = u;
Yd = yd;

Ts = 0.1;
Tend = 10;
Sam_t = 0:Ts:Tend-Ts;
T_length = length(Sam_t);


IDset_ =  struct('MarkovOrder', [8], 'Alpha', [20], 'Beta', [20], 'n', 0.002, 'MinRA', 'era')
[G_ok,H_ok,C_ok,D_ok,F_ok,Sigma, er, M, Ob_CanF] = ...
    OKID_fun_WT02_6(Ud, Yd, IDset_,[]);

x_o(:,1) = pinv(C_ok)*(C*x_initial);
for kk = 1 : T_length
    y_o(:,kk) = C_ok*x_o(:,kk);
    x_o(:,kk+1) = G_ok*x_o(:,kk)  + H_ok*Ud(:,kk) + F_ok*(y_o(:,kk) - Yd(:,kk)) ;
end
error_Y = Yd - y_o;

ii = 0;
for t = 0 : Ts : Tend-Ts
    ii = ii + 1;
    td(:,ii) = ii;
end

% save OKID_esys_without.mat G_ok H_ok C_ok D_ok F_ok
figure
plot(td,Yd,td,y_o,'.-')
legend('y_{d}','y_{o}')