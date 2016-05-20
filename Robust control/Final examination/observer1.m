% <<< Observer >>> %
function xodot = observer1(t, xo, u, Ao, Bo, Co, Ko, y ,ref_c, Kc, Ec)

yo = Co*xo;
% uc = Ec*ref_c - Kc*xo; % Tracker
xodot = Ao*xo + Bo*u - Ko*(yo - y);

end
%  t : time
%  u : actuating signal
%  y : output of plant
% xo : state variables of observer
% Ao : state  matrix of observer
% Bo :  input matrix of observer
% Co : output matrix of observer
% Ko : observer gain
% yo : output of observer
% xodot : time derivative of xo