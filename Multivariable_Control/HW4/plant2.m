function xdott = plant2(t,x,Kc,Ec,A,B,ref,Tf)
persistent k;
if isempty(k)
    k = 0;
end
u = Ec*ref(:,fix(t/Tf)+1)- Kc*x;
if(t<=0.03 && 0)
    subplot(211)
    plot(t,u(1),'bo')
    hold on
    pause(0.1e-7)
    if t >0.02999 
        ax = gca;
%         set(ax,'YLim',[-600 500])
        xlabel('t(s)')
        ylabel('u{_1}(t)')
        set(get(gca,'YLabel'),'Rotation',0);
        legend('u{_1}(t)')
    end
% set(get(gca,'YLabel'),'Position',get(get(gca,'YLabel'),'Position') - [0.0010 0 0])
    subplot(212)
    plot(t,u(2),'bo')
    hold on
    pause(0.1e-7)
    if(t >0.02999)
        ax = gca;
%         set(ax,'YLim',[-600 500])
        xlabel('t(s)')
        ylabel('u{_2}(t)')
        set(get(gca,'YLabel'),'Rotation',0);
        legend('u{_2}(t)')
    end
end
if(t >= k*Tf)
    k = k + 1;
    dlmwrite('u.mat',[t u'],'-append')
end
xdott=A*x+B*( Ec*ref(:,fix(t/Tf)+2)- Kc*x);

end