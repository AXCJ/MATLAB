clc; clear all; close all; 

load lambda;
% lambda = 3*xa*xs/(rhoA*L^3);

scrsz = get(groot,'ScreenSize');

xa = 10;  % actuator located at xa
xs = 2; % sensor located at xs
rhoA = 0.01; % For a slender uniform cylinder with mass per unit length pA
L = 30; % Total length

% simple rigid rotor
% Gr = tf(lambda,[1 0 0]);
% figure('Name', 'Gr' ,'Position',[9 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Left_Top of screen
% bode(Gr);
% title('G_r(s)')
% set(gca, 'XLim' , [0.01 100])
% 
% Phase lead controller
z = 100;
beta = 3.5;
% C1 =20*tf([1 z], [1 z*beta])
% figure('Name', 'C1' ,'Position',[scrsz(3)/3 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Middle_Top of screen
% bode(C1)
% [Gm,Pm] = margin(C1)
% title('Phase Lead Bode diagram','fontsize',30)

% bode(Gr*C1)


% controller 2
% gamma =3.6;
% [b a] = zp2tf(-z, [-beta*z -gamma*beta*z], 1);
% C2 = 5000*tf(b, a)
% figure('Name', 'C2' ,'Position',[scrsz(3)/1.5 scrsz(4)/1.7 scrsz(3)/3 scrsz(4)/3]) % Right_Top of screen
% bode(C2)
% title('C_2(s)')



k = 1e6;
xi = 0.5e-2;
% wn_p1 = [0 0 -10.0+2196j -10.0-2196j];
p1 = [0 0];
p2 = [-10+2196j -10-2196j];
p3 = [-35.4+7081j -35.4-7081j];
p5 = [-73.3+14666j -73.3-14666j];
p6 = [-124.2+24835j -124.2-24835j];
p7 = [-187.2+37444j -187.2-37444j];
p8 = [-261.6+52328j -261.6-52328j];
p9 = [-346.5+ 69307j -346.5-69307j];
% xi_z1 = 0.005;
z1 = [-8.51+1933j -8.51-1933j];
z3 = [-37.7+7307j -37.7-7307j];
z4 = [-62.6+13721j -62.6-13721j];
z5 = [-89.7+21032j -89.7-21032j];
z6 = [-158.5+33656j -158.5-33656j];
z7 = [-248.6+50538j -248.6-50538j];
z9 = [-366.5+71763j -366.5-71763j];

[b a] = zp2tf([z1 z3 z4 z5 z6 z7 z9]', [p1 p2 p3 p5 p6 p7 p8 p9], k);
G = tf(b, a);
figure('Name', 'G_flexible' ,'Position',[9 scrsz(4)/1.7/4 scrsz(3)/3 scrsz(4)/3]) % Left_Bottom of screen
bode(G)
title('G_{flexible}(s)')
hold on
s = tf('s');
% lambda = bode(G,1);
bode(lambda/s^2)
grid












