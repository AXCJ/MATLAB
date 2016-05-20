function dxdt=fun_chaotic_system(t,x)

a=35;b=3;c=28;



dxdt(1)=a*(x(2)-x(1));
dxdt(2)=(c-a)*x(1)-x(1)*x(3)+c*x(2);
dxdt(3)=x(1)*x(2)-b*x(3);

dxdt=dxdt';
end