function d = disturb1( t )

d(: , 1) = sin(2*t);
d(: , 2) = 0.25*(sin(4*pi*t)+cos(2*pi*t)+sin(pi*t)+sin(0.5*pi*t))-1;
% d(: , 2) = cos(2*t);
end

