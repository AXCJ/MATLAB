% Simulation_Scopes
%
% Plot of various oscillograms recorded during simulations
%
%  The names entered for each oscilloscope in the simulation diagrams,
%  in the field "Variable name" which appears under the tab "Data history" after
%  having clicked on the oscilloscope and then on its "Parameters" icon,
%  must matched the names used below, in this program.
%
%  Author: E. Ostertag, 10 February 2010  
%  Last update: 3 March 2011

rep = input('Plot the simulation oscillograms (Y/N)? [N] --> ','s');
if isempty(rep) || upper(rep(1)) == 'N'
  return;
end
%
% Oscillograms of the state variables and (eventually) their estimates
%
if iobs == 0,
  plot_xi = input('Real states (Y/N)? [N] --> ','s');
else
  plot_xi = input('Real and estimated states (Y/N)? [N] --> ','s');
end
if isempty(plot_xi)
  plot_xi = 'N';
else
  plot_xi = upper(plot_xi(1));
end
if plot_xi == 'Y',
  if iobs == 0,
    superimposed_xi = 'N';
  else
    if n >= 3
      superimposed_xi = input('   3 first states, real and estimated curves superimposed (Y/N)? [N] --> ','s');
    elseif n == 2
      superimposed_xi = input('   2 states, real and estimated curves superimposed (Y/N)? [N] --> ','s');
    else
      superimposed_xi = input('   unique state, real and estimated curves superimposed (Y/N)? [N] --> ','s');
    end
  end
  if isempty(superimposed_xi)
    superimposed_xi = 'N';
  else
    superimposed_xi = upper(superimposed_xi(1));
  end
  if superimposed_xi == 'Y',
    figure; set(gcf,'defaultlineLineWidth',1.5);
    if plant == 'C',
      plot(CompareScope.time, CompareScope.signals(1).values(:,1),'r-'); hold on;
      plot(CompareScope.time, CompareScope.signals(1).values(:,2),'b--'); grid on;
    else
      stairs(CompareScope.time, CompareScope.signals(1).values(:,1),'r-'); hold on;
      stairs(CompareScope.time, CompareScope.signals(1).values(:,2),'b--'); grid on;
    end
    legend('real','estimated');
    xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
    title('x1, real and estimated','FontName','Times New Roman','FontSize',14.0);
    if n >= 2,
      figure; set(gcf,'defaultlineLineWidth',1.5);
      if plant == 'C',
        plot(CompareScope.time, CompareScope.signals(2).values(:,1),'r-'); hold on;
        plot(CompareScope.time, CompareScope.signals(2).values(:,2),'b--'); grid on;
      else
        stairs(CompareScope.time, CompareScope.signals(2).values(:,1),'r-'); hold on;
        stairs(CompareScope.time, CompareScope.signals(2).values(:,2),'b--'); grid on;
      end
      xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
      title('x2, real and estimated','FontName','Times New Roman','FontSize',14.0);
      legend('real','estimated');
    end
    if n >= 3
      figure; set(gcf,'defaultlineLineWidth',1.5);
      if plant == 'C',
        plot(StateScope.time, StateScope.signals(3).values,'r-'); hold on;
        plot(EstimateScope.time, EstimateScope.signals(3).values,'b--'); grid on;
      else
        stairs(StateScope.time, StateScope.signals(3).values,'r-'); hold on;
        stairs(EstimateScope.time, EstimateScope.signals(3).values,'b--'); grid on;
      end
      xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
      title('x3, real and estimated','FontName','Times New Roman','FontSize',14.0);
      legend('real','estimated');
    end
  else
    for i = 1:n
      if i > size(StateScope.signals,2)
        disp('*** Not enough scope channels to continue ***'); beep; break
      end
      figure; set(gcf,'defaultlineLineWidth',1.5);
      if plant == 'C',
        plot(StateScope.time, StateScope.signals(i).values,'r-'); grid on;
      else
        stairs(StateScope.time, StateScope.signals(i).values,'r-'); grid on;
      end
      xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
      title(['x' num2str(i) ' real'],'FontName','Times New Roman','FontSize',14.0);
      if iobs > 0
        figure; set(gcf,'defaultlineLineWidth',1.5);
        if plant == 'C',
          plot(EstimateScope.time, EstimateScope.signals(i).values,'b--'); grid on;
        else
          stairs(EstimateScope.time, EstimateScope.signals(i).values,'b--'); grid on;
        end
        xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
        title(['x' num2str(i) ' estimated'],'FontName','Times New Roman','FontSize',14.0);
      end
    end
  end
end
%
% Oscillograms of measured output(s)
%
if q == 1
  plot_yi = input('Output y (Y/N)? [N] --> ','s');
else
  plot_yi = input(['Outputs yi, i=1,' num2str(q) ' (Y/N)? [N] --> '],'s');
end
if isempty(plot_yi)
  plot_yi = 'N';
else
  plot_yi = upper(plot_yi(1));
end
if plot_yi == 'Y',
  for i = 1:q
    if i > size(MeasurementScope.signals,2)
      disp('*** Not enough scope channels to continue ***'); beep; break
    end
    figure; set(gcf,'defaultlineLineWidth',1.5);
    if plant == 'C',
      plot(MeasurementScope.time, MeasurementScope.signals(i).values,'r-'); grid on;
    else
      stairs(MeasurementScope.time, MeasurementScope.signals(i).values,'r-'); grid on;
    end
    xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
    if q == 1,
      title('Output y','FontName','Times New Roman','FontSize',14.0);
    else
      title(['Output y' num2str(i) ],'FontName','Times New Roman','FontSize',14.0);
    end
  end
end
%
% Oscillograms of control signal(s)
%
if p == 1
  plot_ui = input('Control u (Y/N)? [N] --> ','s');
else
  plot_ui = input('Controls ui (i=1,...,p)? [N] --> ','s');
end
if isempty(plot_ui)
  plot_ui = 'N';
else
  plot_ui = upper(plot_ui(1));
end
if plot_ui == 'Y',
  if p == 1
    superimposed_ui = 'N';
  else
	  superimposed_ui = input('   plots superimposed on the same axes (Y/N)? [N] --> ','s');
    if isempty(superimposed_ui)
      superimposed_ui = 'N';
    else
      superimposed_ui = upper(superimposed_ui(1));
    end
  end
	if superimposed_ui == 'Y'
    figure; set(gcf,'defaultlineLineWidth',1.5);
    if plant == 'C',
      plot(ControlScope.time, ControlScope.signals.values(:,1),'r-'); hold on;
    else
      stairs(ControlScope.time, ControlScope.signals.values(:,1),'r-'); hold on;
    end
    if p >= 2,
      if plant == 'C',
        plot(ControlScope.time, ControlScope.signals.values(:,2),'b-'); hold on;
      else
        stairs(ControlScope.time, ControlScope.signals.values(:,2),'b-'); hold on;
      end
    end
    if p >= 3,
      if plant == 'C',
        plot(ControlScope.time, ControlScope.signals.values(:,3),'g-');
      else
        stairs(ControlScope.time, ControlScope.signals.values(:,3),'g-');
      end
    end;
    grid on;
    xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
    if p == 1
      title('Closed-loop control signal u (step response)','FontName','Times New Roman','FontSize',14.0);
    elseif p == 2,
      title('Closed-loop control signals u1 and u2 (step responses)','FontName','Times New Roman','FontSize',14.0);
      legend('u1','u2');
    elseif p == 3,
      title('Closed-loop control signals u1, u2 and u3 (step response)','FontName','Times New Roman','FontSize',14.0);
      legend('u1','u2','u3');
    end
  else
    figure; set(gcf,'defaultlineLineWidth',1.5);
    for j=1:p,
      subplot(p,1,j);
      if plant == 'C',
        plot(ControlScope.time, ControlScope.signals.values(:,j),'r-'); grid on;
      else
        stairs(ControlScope.time, ControlScope.signals.values(:,j),'r-'); grid on;
      end
      if j == 1, 
        title('Closed-loop control signal(s) (step response)','FontName','Times New Roman','FontSize',14.0);
      end
      if p == 1,
        ylabel('u','FontName','Times New Roman','FontSize',12.0);
      else
        ylabel(['u' num2str(j,1) ],'FontName','Times New Roman','FontSize',12.0);
      end
    end
    xlabel(tscale,'FontName','Times New Roman','FontSize',12.0);
  end
end
