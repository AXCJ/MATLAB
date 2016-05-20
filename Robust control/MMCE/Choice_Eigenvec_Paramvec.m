function [V, P, choice, kident] = Choice_Eigenvec_Paramvec(A,B,n,p,Lambda,lambdades)
%
% Acquisiton of eigenvectors and parameter vectors
%
%  [V, P] = Choice_Eigenvec_Paramvec(A,B,n,p,Lambda,lambdades)  compares first the 
%  eigenvalues prescribed for the closed-loop system, transmitted by the input vector 
%  argument lambdades, with the eigenvalues of the open-loop (matrix A),
%  which are transmitted by another input argument, the diagonal matrix Lambda.
%  Two cases are possible:
%
%     * case 1: there is no common eigenvalue to the open-loop and the closed-loop system;
%     closed-loop parameter vectors can then be entered either entirely at random, 
%     or in accordance with the choice of the closed-loop eigenvectors ("targeted" choice);
%     in this case, the subroutine returns "choice = 0" and the vector "kident = 0";
%
%     * case 2: one eigenvalue at least is common; the random choice of the parameter
%     vectors is the only possible one; the subroutine returns then 'choice = 1', and
%     a vector kident which contains as component "i" (corresponding to the
%     i-th closed-loop prescribed eigenvalue) the index "j" of the first eigenvalue of A
%     which is found identical with this value. 
%
%  Author: E. Ostertag, 10 February 2010
%  Last update: 11 August 2010
%

% 1. Looking for closed-loop eigenvalues which are identical to one open-loop eigenvalue
%
choice = 0;
for i=1:n,
  for j=1:size(Lambda,1),
    if lambdades(i) == Lambda(j,j),
      kident(i) = j; choice = 1; break
    else,
      kident(i) = 0;
    end
  end
end
%
% 2. If the closed-loop eigenvalues are all distinct from plant eigenvalues,
%    possibiliy of targeted choice of eigenvectors and parameter vectors
%
if choice == 0,
  choice = input('Choice of parameter vectors: 1 = at random; 2 = targeted; Your choice [1] --> ');
  if isempty(choice),
    choice = 1; disp('1');
  end
end
if choice == 1,
  disp(['Enter n = ' num2str(n) ' parameter vectors of size p = ' num2str(p) ' ===>']);
  for i = 1:n,
    pi = zeros(p+1,1);
    while length(pi) ~= p,
      pi = input(['    Parameter vector p' num2str(i) ' = ']);
      if length(pi) == p,
        if size(pi,1) == 1,
          P(:,i) = pi';
        else,
          P(:,i) = pi;
        end
      else,
        disp(['Enter this vector again (size ' num2str(p) ')']); beep;
      end
    end
  end
  V = zeros(n);
else
  disp('Targeted choice of certain components of the eigenvectors and of the parameter vectors')
  nstrg = num2str(n); pstrg = num2str(p);
  for i=1:n,
    M = [A-lambdades(i)*eye(n) -B]; N = zeros(n,1);
    rM = rank(M); jmax = n+p-rM;
    disp(['Eigenvector/parameter vector associated with lambda = ' ...
          num2str(lambdades(i)) ': ' num2str(jmax) ' unknowns to choose arbitrarily:']);
    disp(['  Eigenvector (' nstrg ' components): unknowns 1 to ' nstrg '; '...
        'parameter vector (' pstrg ' components): unknowns ' num2str(n+1) ' to ' num2str(n+p)]);
    for j=1:jmax
      ni(j) = input('        Number of unknown --> ');
      vi(j) = input('        Arbitrary value chosen --> ');
      N = N - M(:,ni(j)-j+1)*vi(j); M(:, ni(j)-j+1) = []; 
    end
    sol = inv(M)*N;
    [Y,I] = sort(ni);
    for j=1:jmax
      sol = [sol(1:Y(j)-1); vi(I(j)); sol(Y(j):end)];
    end
    V(:,i) = sol(1:n);
    P(:,i) = sol(n+1:end);
  end
end
