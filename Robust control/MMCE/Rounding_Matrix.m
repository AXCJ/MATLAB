function X = Rounding_Matrix(X)
%
% Rounds the elements of a matrix with absolute values smaller than SQRT(EPS) ~= 1.5*10^-8 to 0.
%
%  X = Rounding_Matrix(X)  returns matrix X after rounding to zero its elements
%  with absolute values smaller than SQRT(EPS).
%
%  Author: E. Ostertag, 10 February 2010
%  Last update: 3 March 2011
%

p = size(X,1); q = size(X,2);
for i=1:p,
  for j=1:q,
    if abs(X(i,j)) <= sqrt(eps),
      X(i,j) = 0;
    end
  end
end