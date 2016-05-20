% <<< System State Equation >>> %
function xdot = plant1( t, x, u, A, B)

xdot = A*x + B*u;

end
% t : time
% x : state variables
% u : actuating signal
% A : state matrix
% B : input matrix
% xdot : time derivative of x
