t = 0:0.01:10;
x = 10*t;
y = 5*t;
plot(t,x);
print('-clipboard','-dbitmap');

figure(2);
plot(t,y);
print('-clipboard','-dbitmap');