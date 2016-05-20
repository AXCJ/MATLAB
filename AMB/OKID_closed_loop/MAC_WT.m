function [ M ] = MAC_WT( U, Sigma, V, A, B, C, Alpha, Beta)
%MAC_WT Summary of this function goes here
%   Detailed explanation goes here

n = size(A, 1);
[Ve, Dia] = eig(A);
P_tiu = U*Sigma^(0.5)*Ve;
Cm = C*Ve;
P_bar = [Cm];
for ii = 1:(Alpha-1)
    P_bar = [P_bar;
            Cm*Dia^ii];
end

Bm = inv(Ve)*B;
Q_bar = [Bm];
Q_tiu = [B];
for ii = 1:(Beta-1)
    Q_bar = [Q_bar, (Dia^ii)*Bm];
    Q_tiu = [Q_tiu, A^ii*B];
end
Q_tiu = inv(Ve)*Q_tiu;

for ii = 1:n
    IMAC(ii, :) = abs(Q_tiu(ii, :)*Q_bar(ii, :)')/...
        ( sqrt(Q_tiu(ii, :)*Q_tiu(ii, :)') * sqrt(Q_bar(ii, :)*Q_bar(ii, :)'));
end

for ii = 1:n
    OMAC(ii, :) = abs(P_tiu(:, ii)'*P_bar(:, ii))/...
        ( sqrt(P_tiu(:, ii)'*P_tiu(:, ii)) * sqrt(P_bar(:, ii)'*P_bar(:, ii)));
end

TMAC = IMAC.*OMAC;
eigA = eig(A);
M = [eigA, IMAC, OMAC, TMAC];
end

