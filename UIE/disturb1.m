function d = disturb1( t )

a1 = 2*sqrt(t);
w1 = 0.1*cos(pi*t/5);
a2 = 2*cos(pi*t/5);
w2 = 0.25*t+1;
d(: , 1) = sin(2*t);
d(: , 2) = a1.*sin(w1.*t) + a2.*sin(w2.*t);
% d(: , 2) = cos(2*t);
end

