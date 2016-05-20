%     figure('Renderer','opengl')
figure
    hold on

    h = plot(0,0,'XDataSource','x','YDataSource','y');
%     anno=annotation('arrow',[0 0],[0.1 0.1],'color','r');
    htf = hgtransform('parent',gca);
    fill([0 -0.08 0.15 -0.08]-0.12,[0 0.1 0 -0.1],'r', 'edgecolor','none','parent',htf,'linesmoothing','on', 'zdata',ones(4,1));
    axis([0,2*pi,-3 3])
    axis equal
    %%
    N=200;

    xx=linspace(0,2*pi,N);
    for k = 1:N;
        x=xx(1:k);
        y=sin(x);
        refreshdata(h)
%         M= makehgtform('translate',[x(end) y(end) 0],'zrotate',pi/4*cos(x(end)));
        M= makehgtform('translate',[x(end) y(end) 0],'zrotate', cos(x(end)));
        cosXend = cos(x(end))
        Yend = y(end)
        set(htf,'Matrix',M)
        pause(1e-5)
        drawnow
    end
    


%{
figure
x=0:0.1:10;y=x.^2;h=plot(x,y,'o',x,y);
set(gca,'YTick',[0,10,25,50,80,99],'XTick',[0.5,8,10]);
set(gca,'YTickLabel',{0 10 25 50 'cutoff' 99});

figure
t=0:0.01:2*pi;
y=exp(sin(t));
h=plot(t,y,'YDataSource','y');
for k=1:0.1:20
y=exp(sin(t.*k));
refreshdata(h,'caller');
drawnow;
pause(0.1);
end
%}



figure
Tf = 0.02;
sigma = sqrt(35/12); % standard deviation
mean = 3.5; % mean
x = 0:Tf:10;
% bell-shaped curve
fx = @(x) (sigma*sqrt(2*pi))\exp(-((x-mean).^2)/2*sigma^2);
f=fx(x);
plot(x,f)
h = gca;
xtick = [ mean mean+sigma mean-sigma];
set(h, 'XTick', sort(xtick))
grid

L = linspace(-5,5,100);
plot(sin(L)); hold on
for i = 1:length(L)
    j(i) = sat(sin(L(i)));
end
plot(j)

