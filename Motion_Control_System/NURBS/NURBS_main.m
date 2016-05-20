%{   
NURBS

k : order ¦¸¼Æ
p = k-1 : degree ¶¥¼Æ
W : Weights 
number of control points : n+1
number of knots : n+k+1 = r+1
%}

clc;clear;close all;
Linewidth = 3; % Line width in figure
FontSize = 25; % Font size in figure
if 1
    isNFigureVisible = 'on';
else
    isNFigureVisible = 'off';
end
if 1
    isRFigureVisible = 'on';
else
    isRFigureVisible = 'off';
end

p = 4;
K = p+1;
u = [0 0 0 0 0 0.5 1 1 1 1 1]; % knots (Target)
W = [1 1 1 1.5 1 1]; % weights (Target)
% A = [0.0 0.0 0.0 1.0 1.0 1.0]; % control points (Target)
r = length(u) - 1;
W = [W zeros(1 , r - length(W))];
n = r - K;
L = linspace(u(1),u(end),r*100); % number of sampling data
Nfig_name = char('N_{j,1}' , 'N_{j,2}' , 'N_{j,3}' , 'N_{j,4}' , 'N_{j,5}' , 'N_{j,6}');
Rfig_name = char('R_{j,1}' , 'R_{j,2}' , 'R_{j,3}' , 'R_{j,4}' , 'R_{j,5}' , 'R_{j,6}');
Cfig_name = char('Position' , 'Velocity' , 'Acceleration');

% ~~~ N First-order
kk = 1; % order index
figure('Name',Nfig_name(1,:) , 'visible' , isNFigureVisible) % figure 1
hold on
for i = 1 : r % u(0)~u(r-1)
    for j = 1 : length(L)
        if L(j) == u(end)
            N(i , j , kk) = N(i , j-1 , kk);
        else
            if u(i) <= L(j) && L(j) < u(i+1)
                N(i , j , kk) = 1;
            else
                N(i , j , kk) = 0;
            end
        end
    end
    plot(L,N(i , : , kk),'LineWidth' , Linewidth)
end
xlabel(Nfig_name(1,:) , 'FontSize' , FontSize , 'fontweight' , 'b');
hold off
% ~~~ First-order ~~~End


% ~~~ N start at second-order ~~~
for k = 1 : K-1
    kk = kk+1; % start at kk = 2
    figure('Name' , Nfig_name(kk,:) , 'visible' , isNFigureVisible)
    hold on
    for i = 1 : r-k % perform all of N
        for j = 1 : length(L) % perform each N
            N(i , j , kk) = Nik_fun(L(j) , u(i) , u(i+kk-1) , N(i , j , kk-1) , ...
                                  u(i+kk) , u(i+1) , N(i+1 , j , kk-1));
        end
        plot(L , N(i , : , kk) , 'LineWidth' , Linewidth)
    end
    xlabel(Nfig_name(kk,:) , 'FontSize' , FontSize , 'fontweight' , 'b');
    hold off
end
% ~~~ N ~~~ End

% ~~~ N' & N'' 
for kk = 2 : K
    for i = 1 : r-(kk-1)
        Ndot(i , : , kk) = Nik_md_fun(u(i) , u(i+kk-1) , N(i , : , kk-1) , ...
                                                      u(i+kk) , u(i+1) , N(i+1 , : , kk-1) , kk);
    end
    Ndot(i+1 , : , kk) = 0;
    for i = 1 : r-(kk-1)
        Ndot2(i , : , kk) = Nik_md_fun(u(i) , u(i+kk-1) , Ndot(i , : , kk-1) , ...
                                                      u(i+kk) , u(i+1) , Ndot(i+1 , : , kk-1) , kk);
    end
end
% ~~~ N' & N'' ~~~ End

% ~~~ R & R' & R''~~~ 
kk = 0;
checkR(1 , : , :) = 0;
for k = 1 : K
    kk = kk + 1;
    ii = 0; 
%     if(kk == 1)
        figure('Name',Rfig_name(kk,:) , 'visible' , isRFigureVisible)
%     end
    hold on
    
    WN(: , : , kk) = W * N(: , : , kk); % perform WN
    WNdot(: , : , kk) = W(1:size(Ndot,1)) * Ndot(: , : , kk); % perform WN'
    WNdot2(: , : , kk) = W(1:size(Ndot2,1)) * Ndot2(: , : , kk); % perform WN''
    
    for i = 0 : n
        ii = ii + 1;
        R(ii, : , kk) = (W(ii) * N(ii , : , kk)) ./ WN(: , : , kk); % ~~~ perform R
        plot(L,R(ii , : , kk), 'LineWidth' , Linewidth)
        % ~~~ perform R'
        Rdot(ii, : , kk) = W(ii) * Ndot(ii , : , kk) ./ WN(: , : , kk) - ...
                                    W(ii) * N(ii , : , kk) .* WNdot(: , : , kk) ./ WN(: , : , kk) .^ 2;
        % ~~~ perform R''
        Rdot2(ii, : , kk) = W(ii) * Ndot2(ii , : , kk) ./ WN(: , : , kk) - ...
                                     (2 * W(ii) * Ndot(ii , : , kk) .* WNdot(: , : , kk) + ...
                                     W(ii) * N(ii , : , kk) .* WNdot2(: , : , kk)) ./ (WN(: , : , kk) .^ 2) + ...
                                     (2 * W(ii) * N(ii , : , kk) .* (WNdot(: , : , kk) .^ 2)) ./ (WN(: , : , kk) .^ 3);
    end
    xlabel(Rfig_name(kk,:) , 'FontSize' , FontSize , 'fontweight' , 'b');
    hold off
end
% ~~~ R & R' & R'' ~~~ End
% ~~~ check R
sumR = sum(R , 1);
sumRdot = sum(Rdot , 1);
sumRdot2 = sum(Rdot2 , 1);

if(find(sumR < 0.999999 | sumR >1.00001))
    fprintf('-------Error-------R is fault\n')
end
if(find(sumRdot < -0.0000001 | sumRdot >0.00001))
    fprintf('-------Error-------Rdot is fault\n')
end
if(find(sumRdot2 < -0.0000001 | sumRdot2 >0.00001))
    fprintf('-------Error-------Rdot2 is fault\n')
end
% ~~~ check R ~~~ End 

% Calculate control point with motion constraints
% cond_R*A = cond_C => A = inv(cond_R)*cond_C
A = (inv([R( : , 1 , end) Rdot( : , 1 , end) Rdot2( : , 1 , end) ...
                R( : , end , end) Rdot( : , end , end) Rdot2( : , end , end)]')*[0 0 0 1 0 0]')';
% ~~~ Perform control point with motion constraints ~~~ End

% ~~~ C ~~~
for kk = 1 : K
    sumC = 0; sumCdot = 0; sumCdot2 = 0; ii = 0;
    for i = 0 : n
        ii = ii + 1;
        
        subC(ii , : , kk) = A(: , ii) * R(ii , : , kk);
        sumC= sumC + subC(ii , : , kk);

        subCdot(ii , : , kk) = A(: , ii) * Rdot(ii , : , kk);
        sumCdot= sumCdot + subCdot(ii , : , kk);
        
        subCdot2(ii , : , kk) = A(: , ii) * Rdot2(ii , : , kk);
        sumCdot2= sumCdot2 + subCdot2(ii , : , kk);
    end
    C(: , : , kk) = sumC;
    Cdot(: , : , kk) = sumCdot;
    Cdot2(: , : , kk) = sumCdot2;
end
% for plot
Cprime(1 , : ) = C(: , : , kk);
Cprime(2 , : ) = Cdot(: , : , kk);
Cprime(3 , : ) = Cdot2(: , : , kk);
% ~~~ C ~~~ End

for i = 1 : 3
    figure('Name',Cfig_name(i,:) , 'visible' , isNFigureVisible)
    plot(L,Cprime(i,:),'LineWidth',Linewidth)
    xlabel(Cfig_name(i,:) , 'FontSize' , FontSize , 'fontweight' , 'b');
end

h1 = figure('Name' ,'C & C'' & C''''' , 'Position' , [20 558 560 420]);
plot(L , C(: , : , 5)*10 ,'r', L , Cdot(: , : , 5)*5  ,'b' , L , Cdot2(: , : , 5) ,'k', 'LineWidth' , Linewidth)
ax = gca;
set(ax,'box','off' ,'XTick',[0 1],'YTick',[],'YLim',[-6 10],'FontSize' , FontSize)
text(L(length(L)*0.8),C(: , length(L)*0.9 , 5)*9,'\color{red} C(u)',...
     'HorizontalAlignment','left',...
     'FontSize',FontSize)
text(L(length(L)*0.7),Cdot(: , length(L)*0.7 , 5)*5+0.4,'\color{Blue} C^{(1)} (u)',...
     'HorizontalAlignment','left',...
     'FontSize',FontSize)
text(L(length(L)*0.6),Cdot2(: , length(L)*0.7 , 5)+2.2,'\color{black} C^{(2)} (u)',...
     'HorizontalAlignment','left',...
     'FontSize',FontSize)

