%*************************************************************************
%  LQ_robustness_continuous:
%    Study of the robustness of the continuous-time LQ control law
%*************************************************************************
%
%  Author: E. Ostertag, 10 February 2010
%  Last update: 11 August 2010
%
echo on
%
% 1) Plot circle which is "forbidden" to the Nyquist diagrams
%
echo off
figure; theta = 0:2*pi/360:2*pi;
set(gca,'DefaultLineLineWidth',1.5);
axis([-5 1 -2 2]);
axis equal; axis manual; hold on; grid on;
plot(-1 + cos(theta),sin(theta),'r'); % Plot "forbidden" circle
echo on
%
% 2) Plot on the same axes the Nyquist diagram of the compensated open loop; 
%    suitable angular frequency range: 0.1 rad/s to 100 rad/s
%
echo off
Topen_ss = ss(A,B,L,0);      % Compensated open loop
nbpoints = 400; decmin = -1; decmax = 2;
omega = logspace(decmin,decmax,nbpoints);
[Re,Im] = nyquist(Topen_ss,omega);    % Compute Nyquist diagram
Re_Topen = zeros(1,nbpoints); Im_Topen = zeros(1,nbpoints);
Re_Topen(:) = Re(1,1,:); Im_Topen(:) = Im(1,1,:);
plot(Re_Topen,Im_Topen,'b'); plot(Re_Topen,-Im_Topen,'m');  % Plot Nyquist diagram
legend('"forbidden" circle','compensated open-loop Nyquist diagram, omega > 0','compensated open-loop Nyquist diagram, omega < 0');
