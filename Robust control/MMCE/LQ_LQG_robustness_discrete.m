%*************************************************************************
%  LQ_LQG_robustness_discrete:
%    Study of the robustness of the LQ and the LQG control law, in the
%    discrete-time case
%*************************************************************************
%
%  Author: E. Ostertag, 10 February 2010
%  Last update: 17 August 2010
%
if isempty(KEST)
  echo on
% 1) Study of the LQ control robustness (without filter)
%
% a) Plot the unit circle and the circle "forbidden" to the
%    Nyquist diagram
%
  echo off
  figure; theta = 0:2*pi/360:2*pi;
  set(gca,'DefaultLineLineWidth',1.5);
  axis([-2 1 -1 1]);
  axis equal; axis manual;  hold on;
  plot(cos(theta),sin(theta),'g');      % Plot unit circle
  rho = 1/sqrt(1+(1/R)*Gamma'*P*Gamma);
  plot(-1 + rho*cos(theta),rho*sin(theta),'r'); % Plot "forbidden" circle
  echo on
%
% b) Plot on the same axes the compensated open loop diagram,
%    without filter in the loop;
%    suitable angular frequency range: 0.01 rad/s to pi/Ts
%
  echo off
  Topen_ss = ss(Phi,Gamma,L,D,Ts);      % Compensated open loop
  nbpoints = 400; decmin = -2; decmax = log10(pi/Ts);
  omega = logspace(decmin,decmax,nbpoints);
  [Re,Im] = nyquist(Topen_ss,omega);    % Compute Nyquist diagram
  Re_Topen = zeros(1,nbpoints); Im_Topen = zeros(1,nbpoints);
  Re_Topen(:) = Re(1,1,:); Im_Topen(:) = Im(1,1,:);
  plot(Re_Topen,Im_Topen,'b');   grid on;  % Plot Nyquist diagram
  legend('unit circle','"forbidden" circle','compensated open-loop Nyquist diagram');
end
%
if ~isempty(KEST)
  echo on
% 2) Study of the LQG control robustness (filter in the loop)
%
% a) Plot unit circle and "forbidden" circle
%
  echo off
  figure; theta = 0:2*pi/360:2*pi;
  set(gca,'DefaultLineLineWidth',1.5);
  axis([-2 1 -1 1]);
  axis equal; axis manual; grid on; hold on;
  plot(cos(theta),sin(theta),'g');
  rho = 1/sqrt(1+(1/R)*Gamma'*P*Gamma);
  plot(-1 + rho*cos(theta),rho*sin(theta),'r');
  echo on
%
% b) Plot on the same axes the Nyquist diagram of the compensated
%    discrete-time open loop;
%    suitable angular frequency range: 0.01 rad/s to pi/Ts
%
  echo off
  RLQG = lqgreg(KEST,L,'current');
  PT326_ext = ss(Phi, Gamma, Cmm, D, Ts)
  Topen_ss = series(PT326_ext, RLQG);
  nbpoints = 400; decmin = -2; decmax = log10(pi/Ts);
  omega = logspace(decmin,decmax,nbpoints);
  [Re,Im] = nyquist(Topen_ss,omega);
  Re_Topen = zeros(1,nbpoints); Im_Topen = zeros(1,nbpoints);
  Re_Topen(:) = Re(1,1,:); Im_Topen(:) = Im(1,1,:);
  plot(Re_Topen,Im_Topen,'b');
  legend('unit circle','"forbidden" circle','compensated open-loop Nyquist diagram');
  echo on
end
