% <<< System State Equation >>> %
function xdot = plant1( t, x, u, A, B)

xdot = A*x + B*u;

end

