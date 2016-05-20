function [A, B, C, D, G, Sigma, er, M, Ob_CanF] = OKID_fun_WT02_6(u, y, IDdesign_set, dd)
% IDdesign_set = struct('MarkovOrder', [1], 'Alpha', [40], 'Beta', [20], 'n', 0.03);
% m: output number
% r: input number
% p: Markov parameters order
% l: Sample data points
Alpha = IDdesign_set.Alpha;
Beta = IDdesign_set.Beta;
p = IDdesign_set.MarkovOrder;
n = IDdesign_set.n;

[m,l] = size(y);
r = size(u,1);

N = Alpha + Beta;

%===== Get Markov parameters 
upsilon = [ u;
            y];
index_p = p + 1;    
y_bar = [y(:, index_p:l)];
V_bar = [u(:, index_p:l)];
% if nargin == 4
%     V_bar = [zeros(size(u(:, index_p:l)))];
% end
if nargin == 4
    V_bar = [];
    D = zeros(size(y, 1), size(u, 1));
end

for ii = 1:p
    V_bar = [V_bar;
             upsilon(:, (index_p-ii):(l-ii) )];
end
Y_bar = y_bar*pinv(V_bar);
if nargin ~= 4
    D =  Y_bar(:, 1:r); % D : ID systeim => D
    Y_bar(:, 1:r) = [];
end

for ii = 1:p
    Yk_bar(:, :, ii) = Y_bar(:, ((ii-1)*(r+m)+1):(ii*(r+m)) );
end

for jj = 1:length(Yk_bar(1, 1, :))
    Yk1_bar(:, :, jj) = Yk_bar(:, 1:r, jj);
    Yk2_bar(:, :, jj) = -Yk_bar(:, ((r+1):(r+m)), jj);
end

Pk(:, :, 1) = [ (Yk1_bar(:, :, 1)-Yk2_bar(:, :, 1)*D),  Yk2_bar(:, :, 1)];
Yk(:, :, 1) = Pk(:, 1:r, 1);
Yok(:, :, 1) = Pk(:, (r+1):(r+m), 1);

for kk = 2:N
    if (kk <= p)
        sumTempPk = zeros(size(Pk(:, :, 1) ));
        for ii = 1:(kk-1)
            sumTempPk = sumTempPk + Yk2_bar(:, :, ii)*Pk(:, :, kk-ii);
        end
        Pk(:, :, kk) = [ (Yk1_bar(:, :, kk)-Yk2_bar(:, :, kk)*D),  Yk2_bar(:, :, kk)] - sumTempPk;
    else
        sumTempPk = zeros(size(Pk(:, :, 1) ));
        for ii = 1:p
            sumTempPk = sumTempPk + Yk2_bar(:, :, ii)*Pk(:, :, kk-ii);
        end
        Pk(:, :, kk) = -sumTempPk;
    end
    
    Yk(:, :, kk) = Pk(:, 1:r, kk);
    Yok(:, :, kk) = Pk(:, (r+1):(r+m), kk);
end
%===== Get Markov parameters ===end
%===== Hankel matrix
H =  [];
for ii = 1:Alpha
    Htemp = [];
    for jj = 1:(Beta+1)
        Htemp = [Htemp Pk(:, :, ii+jj-1)];
    end
    H = [H;
        Htemp];
end
H0 = H;
H0(:, (end-(r+m)+1):end) = [];
H1 = H;
H1(:, 1:(r+m) ) = [];
        %===== Hankel matrix ===end
switch (IDdesign_set.MinRA)
    case 'eradc'
        [A, Btemp, C, er, U, Sigma, V] = eradc_WT( H1, H0, (r+m), m, n); 
    otherwise
        [A, Btemp, C, er, U, Sigma, V] = era_WT( H1, H0, (r+m), m, n);
end
B = Btemp(:, 1:r); % B : ID systeim => B
G = Btemp(:, (r+1):(r+m) ); % G : ID systeim => G

[ M ] = MAC_WT( U, Sigma, V, A, B, C, Alpha, Beta);

%====
Ob_CanF = struct('A', [], 'B', [], 'C', [], 'D', [], 'G', []);
clear temp;
Atemp = []; Btemp = []; Ctemp = []; Gtemp = [];
for ii = 1:p-1
    Atemp = [Atemp, zeros(m)];
end
Atemp = [Atemp, -Yk2_bar(:, :, p)];
Btemp = [Yk1_bar(:, :, p)-Yk2_bar(:, :, p)*D];
Gtemp = [Yk2_bar(:, :, p)];
Ctemp = [eye(m)];
for ii = 1:p-1
    temp = [];
    for jj = 1:p-1
        if (ii == jj)
            temp = [temp, eye(m)];
        else
            temp = [temp, zeros(m)];
        end
    end
    temp = [temp, -Yk2_bar(:, :, (p-ii))];
    Atemp = [Atemp;
            temp];
    Btemp = [Btemp;
             Yk1_bar(:, :, (p-ii))-Yk2_bar(:, :, (p-ii))*D];
    Gtemp = [Gtemp;
             Yk2_bar(:, :, (p-ii))];
    Ctemp = [zeros(m), Ctemp];
end
Ob_CanF.A = Atemp;
Ob_CanF.B = Btemp;
Ob_CanF.C = Ctemp;
Ob_CanF.D = D;
Ob_CanF.G = Gtemp;

