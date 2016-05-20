% Responses:
%  Script file to plot various free or step responses, of real and/or estimated
%  states, outputs and controls. 
%
%  Author: E. Ostertag, 10 February 2010
%  Last update: 3 March 2011
%

n=size(A_CL,1); p = size(B,2); q =size(Cm,1);
%
% Plot of input-output step responses of the closed-loop system
%
rep = input('Plot input-output closed-loop system step responses, with x0 = 0 (Y/N)? [N] --> ','s');
if isempty(rep)
  rep = 'N'; disp('N');
else
  rep = upper(rep(1));
end
if rep == 'Y',
%
% Step functions applied individually to each reference input, one at a time at t=0
%
  for i=1:p;
    sys_CLy{i} = ss(A_CL,B_CL(:,i),C_CL,D_CL(:,i),T);
    yr_sim = ones(length(t),1);
    y12 = lsim(sys_CLy{i},yr_sim*yr(i),t);  % Inputs applied individually
    for j=1:q,
      figure; plot(t,y12(:,j)); grid on;
      xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
      title(['Closed-loop I/O step response: from yr' num2str(i) ' to y' num2str(j)],'FontName','Times New Roman','FontSize',14.0);
    end
  end
%
% If more than 1 input: Step functions applied simultaneously to all reference inputs, at times
% specified by "start_reference"
%
  if p == 1
    rep = input(['Single input system : \n'...
      '   plot closed-loop control signal in response to reference step (x0 = 0) (Y/N)? [N] --> '],'s');
    if isempty(rep)
      rep = 'N'; disp('N');
    else
      rep = upper(rep(1));
    end
    if rep == 'Y',
      figure
      sys_CL_u = ss(A_CL,B_CL,-L,0,T);
      yr_sim = ones(length(t),q);
      U = lsim(sys_CL_u,yr_sim,t);
      U=U+ yr_sim*(M*Sc)';
      plot(t,U); grid on;
      title('Control signal of closed-loop system (step response)','FontName','Times New Roman','FontSize',14.0);
      ylabel('u','FontName','Times New Roman','FontSize',12.0);
      xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
    end
  else
	  yr_sim = zeros(length(t),q);
    for i=1:q
      t0 = find(t >= start_reference(i));  % Create steps at times "start_reference"
      yr_sim(t0,i) = yr(i);
    end
    ytot = lsim(sys_CL,yr_sim,t);  % Inputs applied simultaneously, at times "start_reference"
    figure;
    for j=1:q;
      subplot(q,1,j); plot(t,ytot(:,j)); grid on;
      if j == 1, 
        title('Closed-loop I/O step response: from yr to ...','FontName','Times New Roman','FontSize',14.0);
      end
      if q == 1,
        ylabel('y','FontName','Times New Roman','FontSize',12.0);
      else
        ylabel(['y' num2str(j,1) ],'FontName','Times New Roman','FontSize',12.0);
      end
    end
    xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
%
% If p > 1: Plot of closed-loop control responses to reference steps
% applied at times specified by "start-reference", as above
%
    rep = input('Closed-loop control signals in response to reference steps (x0 = 0) (Y/N)? [N] --> ','s');
    if isempty(rep)
      rep = 'N'; disp('N');
    else
      rep = upper(rep(1));
    end
    if rep == 'Y',
      figure
      sys_CL_u = ss(A_CL,B_CL,-L,zeros(p,p),T);
      U = lsim(sys_CL_u,yr_sim,t) + yr_sim*(M*Sc)';
%        N.B.: lsim produces a matrix with LENGTH(T) rows and as many
%        columns as outputs in sys_CL; this is why the second contribution
%        to the control sihnal, i.e. M*Sc*yr_sim', has been transposed
      for j=1:p;
       subplot(p,1,j); plot(t,U(:,j)); grid on;
        if j == 1, 
          title('Control signal(s) of closed-loop system (step response)','FontName','Times New Roman','FontSize',14.0);
        end
        ylabel(['u' num2str(j,1) ],'FontName','Times New Roman','FontSize',12.0);
      end
      xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
    end
  end
end
disp('Non-vanishing initial state:'); x0 = xx0
%
% Plot of output free responses of the open-loop system
%
if integ == 1;
  rep = 'N';
else
  rep = input('Free responses of the open-loop system outputs (u = 0; x0 = as above) (Y/N)? [N] --> ','s');
  if isempty(rep)
    rep = 'N'; disp('N');
  else
    rep = upper(rep(1));
  end
end
if rep == 'Y',
  figure;
  sys_OL = ss(A,B,Cm,D,T);
  [Y,t,x] = initial(sys_OL,x0,t);
  for j=1:q,
    subplot(q,1,j); plot(t,Y(:,j)); grid on;
	  if j == 1,
      title('Free response of open-loop system output(s)','FontName','Times New Roman','FontSize',14.0);
    end
    if q == 1,
      ylabel('y','FontName','Times New Roman','FontSize',12.0);
    else
      ylabel(['y' num2str(j,1) ],'FontName','Times New Roman','FontSize',12.0);
    end
  end
  xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
end
%
% Plot of output free responses of the closed-loop system
%
rep = input('Free responses of the closed-loop system outputs (yr = 0; x0 = as above) (Y/N)? [N] --> ','s');
if isempty(rep)
  rep = 'N'; disp('N');
else
  rep = upper(rep(1));
end
if rep == 'Y',
  figure;
  [Y,t,x] = initial(sys_CL,[x0' zeros(1,qi)],t);
  for j = 1:q,
    subplot(q,1,j); plot(t,Y(:,j)); grid on;
	  if j == 1, 
      title('Free response of closed-loop system output(s)','FontName','Times New Roman','FontSize',14.0);
    end
    if q == 1,
      ylabel('y','FontName','Times New Roman','FontSize',12.0);
    else
      ylabel(['y' num2str(j,1) ],'FontName','Times New Roman','FontSize',12.0);
    end
  end
  xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
end
%
% Plot of the control signals of the closed-loop system, in free response from x0
%
rep = input('Control input(s) of closed-loop system (free response from x0 above) (Y/N)? [N] --> ','s');
if isempty(rep)
  rep = 'N'; disp('N');
else
  rep = upper(rep(1));
end
if rep == 'Y',
  figure
  sys_CL_u = ss(A_CL,B_CL,-L,zeros(p,p),T);
  [U,t,x] = initial(sys_CL_u,[x0' zeros(1,qi)],t);
  for j = 1:p,
    subplot(p,1,j); plot(t,U(:,j)); grid on;
    if j == 1, 
      title('Control signal(s) of closed-loop system (free response from x0)','FontName','Times New Roman','FontSize',14.0);
    end
    if p == 1,
      ylabel('u','FontName','Times New Roman','FontSize',12.0);
    else
      ylabel(['u' num2str(j,1) ],'FontName','Times New Roman','FontSize',12.0);
    end
  end
  xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
end
