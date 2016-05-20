% <<< Observer >>> %
function xodot = observer1(t, xo, u, Ao, Bo, Co, Ko, y)

yo = Co*xo;
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