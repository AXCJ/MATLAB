function [A, B, C, er, U, Sigma, V] = era_WT( H1, H0, r, m, n)
%ERA_WT Summary of this function goes here
%   Detailed explanation goes here
%   r:input
%   m:output
%   n:state number or error
%
[Utemp, S, Vtemp] = svd(H0);
sum_S = 0; sum_n = 0;
for ii = 1:min(size(S))
    sum_S = sum_S + S(ii, ii);
end

er = 1; 

U = Utemp(:, 1:n);
Sigma = S(1:n, 1:n);
V = Vtemp(:, 1:n);

A = Sigma^(-0.5)*U'*H1*V*Sigma^(-0.5); % A : ID systeim => A
Btemp = Sigma^(0.5)*V';
Ctemp = U*Sigma^(0.5);
B = Btemp(:, 1:r);
C = Ctemp(1:m, :);

end

