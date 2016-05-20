function [L, M] = Falb_Wolovich_or_Roppenecker(A,B,C,Lambda,Vp,Ts)
%
% Design of a decoupling multivariable control law, by the Falb-Wolovich method,
% or of a complete modal control (Roppenecker): applies only if p = q.
%
%  [L, M] = Falb_Wolovich_or_Roppenecker(A,B,C,Lambda,Vp,Ts) computes, according to
%  one of the two following methods:
%      * Falb-Wolovich method (choice F),
%      * complete modal control method, or Roppenecker's method (choice R),
%  the state-feedback matrix L which yields a closed-loop system in the form of
%  SISO subsystems in parallel. The open-loop eigenvalues which can be shifted to
%  new values prescribed for the closed-loop system depend on the order difference
%  of the initial plant, and are mentioned to the user for him to choose these new values.
%
%  An input value Ts > 0 indicates that these values should be converted to
%  their discrete equivalent at the period Ts, by exp(value*Ts).
%
%  Vp is the matrix whose columns are the plant eigenvectors, for the case one of the
%  desired eigenvalues is identical to one of the plant eigenvalues (non shifted eigenvalue).
%
%  The returned matrix M is the feedforward matrix which provides a unit steady-state
%  gain to the closed-loop system and which is also computed by this algorithm.
%
%  Called subroutines: Get_Eigenvalue, Rounding_Matrix.
%
%  Author: E. Ostertag, 10 February 2010
%  Last updated: 13 February 2011
%

n = size(A,1); p = size(B,2); q = size(C,1); D = zeros(q,p);
if p ~= q,
  disp('*** The decoupling method applies only if p=q ***'); beep;
  L=0; M=0;
  return
end
%
% 1.- Determining the order difference delta(i) of the various outputs yi
%
for i=1:q,
  di = 0;
  Mtest = C(i,:)*A^di*B;
  while C(i,:)*A^di*B == 0,
    di = di+1;
  end
  delta(i)=di+1;
end
%
% 2.- Global order difference "sumdelta" of the multivariable system
%
Delta = 0;
for i=1:q,
  Delta = Delta+delta(i);
end
%
% 3.- Decoupling only possible if the matrix D* computed below is invertible
%
for i=1:q,
  Dstar(i,:) = C(i,:)*A^(delta(i)-1)*B;
end
%
% 4.- Choice of method
%
  disp('Choice of method:');
  disp('                 F = Falb-Wolovich method (default)');
  disp('                 R = Roppenecker''s method (decoupling by complete modal control)');
  method = input('Your choice --> ','s');
  if isempty(method)
    method = 'F'; disp('F');
  else
    method = upper(method(1));
  end
if method == 'F' && det(Dstar) == 0,
  disp(['D* is not invertible ===> '...
        'decoupling not possible via Falb-Wolovich method']); beep;
  ropp = input('Try Roppenecker''s method (Y/N) ? [Y] --> ','s');
  if upper(ropp(1)) == 'Y'
    method = 'R';
  else
    L=0; M=0;
    return
  end
end

sys_OL = ss(A, B, C, D, Ts);
G_OL = tf(sys_OL);
kP = 1;
if Delta == n,
  disp('Delta = n, all eigenvalues can therefore be shifted arbitrarily');
else
  disp(['Delta < n: thus only part of the eigenvalues can be shifted arbitrarily']);
end
  
if method == 'F',                      % Method of Falb-Wolovich
  for i=1:q,
    disp(' ')
    if delta(i) == 1,
        disp(['Shifting unique eigenvalue corresponding to subsystem # ' num2str(i) ' :']);
    else
        disp(['Shifting ' num2str(delta(i)) ' eigenvalues, corresponding to subsystem # ' num2str(i) ' :']);
    end
    disp(' ')
    lambdades = Get_Eigenvalue(delta(i),Ts);
    linep{i} = real(poly(lambdades));
    if Ts > 0,
      knum(i) = polyval(linep{i},1);
    else
      knum(i) = polyval(linep{i},0);
    end
    P(i,:) = C(i,:)*polyvalm(linep{i},A);
  end
  L_decoupling = inv(Dstar)*P;
else                            % Method of Roppenecker (complete modal control)
  Z_OL = zero(sys_OL);
  nz_OL = length(Z_OL);
  if nz_OL == n - Delta,
    Z_OL = Rounding_Matrix(Z_OL);
    i = 1; nz_OL_inst = 0; Z_OL_inst = [];
    while i <= nz_OL
      if real(Z_OL(i)) >= 0,
        Z_OL_inst = [Z_OL_inst; Z_OL(i)];
        Z_OL(i) = [];
        nz_OL_inst = nz_OL_inst + 1;
        nz_OL = nz_OL - 1;
      else
        i = i+1;
      end
    end
%
% Automatic placement of (nz_OLO <= n-Delta) eigenvalues on the "stable"
% zeros of the plant
%
    disp('*')
    disp(['*** Automatic placement of ' num2str(nz_OL) ' eigenvalue(s) on the "stable" zero(s) of the plant ***'])
    disp('*')
    i = 1;
    while i<= nz_OL
      disp(['Compensation of "stable" zero # ' num2str(i) ' = ' num2str(Z_OL(i)) ]);
      Rosen = [A-Z_OL(i)*eye(n) B; C D];
      if isreal(Z_OL(i))
        nsp = null(Rosen,'r'); nv = size(nsp,2);
        if isempty(nsp)
          nsp = null(Rosen); nv = size(nsp,2);
        end
        nsp = Rounding_Matrix(nsp);
        Xz = nsp(1:n,1:nv), Uz = nsp(n+1:n+p,1:nv),
        V(:,kP:kP+nv-1) = Xz;  P(:,kP:kP+nv-1) = -Uz;
        kP = kP+nv;
        i = i+1;
      else
        nsp = null(Rosen); nv = size(nsp,2);
        nsp = Rounding_Matrix(nsp);
        Xz = nsp(1:n,1:nv), Uz = nsp(n+1:n+p,1:nv),
        V(:,kP:kP+nv-1+1) = [real(Xz) imag(Xz)];  P(:,kP:kP+nv-1+1) = -[real(Uz) imag(Uz)];
        kP = kP+nv+1;
        i = i+2;
      end
      if kP > n,
        break;
      end
    end
%
% Placing arbitrarily Delta eigenvalues
%
    disp('*')
    disp('*** Placement of Delta eigenvalues ***')
    disp('*')
    for i=1:q,
      disp(' ')
      if delta(i) == 1,
        disp(['Enter unique eigenvalue corresponding to subsystem # ' num2str(i) ':']);
      else
        disp(['Enter ' num2str(delta(i)) ' eigenvalues, corresponding to subsystem # ' num2str(i) ':']);
      end
      disp(' ')
      lambdades = Get_Eigenvalue(delta(i),Ts);
      linep{i} = real(poly(lambdades));
      if Ts > 0,
        knum(i) = polyval(linep{i},1);
      else
        knum(i) = polyval(linep{i},0);
      end
      e_vect =[zeros(1,i-1) 1 zeros(1,q-i)]';
      for k=1:delta(i),
% Looking for closed-loop eigenvalues which are identical to one open-loop
% eigenvalue (these eigenvalues are not shifted)
        for j=1:size(Lambda,1),
          if lambdades(k) == Lambda(j,j),
            kident(k) = j; break
          else,
            kident(k) = 0;
          end
        end
        if kident(k) == 0
          G_OL_lambdades = evalfr(G_OL,lambdades(k));
          P(:,kP) = -inv(G_OL_lambdades)*e_vect;
          V(:,kP) = inv(A-lambdades(k)*eye(n))*B*P(:,kP);
        else
          disp('*** Warning : chosen eigenvalue = an open-loop eigenvalue ***')
          disp('*** Corresponding parameter vector = 0, and open-loop eigenvector kept ***'); beep
          P(:,kP) = 0; V(:,kP) = Vp(:,kident(k))
        end
        kP = kP+1;
      end
    end
%
% Arbitrary placement of eigenvalues, one for each plant "RHP" zero
%
    if nz_OL_inst > 0
      disp('*')
      disp(['*** Arbitrary placement of ' num2str(nz_OL_inst) ' additional eigenvalue(s), corresponding to the number'])
      disp('     of plant "RHP" zeros')
      disp('*')
      lambdades = Get_Eigenvalue(nz_OL_inst,Ts);
      i=1; k=1;
      while k <= nz_OL_inst
% % %       disp(['"RHP" zero # ' num2str(i) ' = ' num2str(Z_OL_inst(i)) ]); 
        e_vect =[zeros(1,i-1) 1 zeros(1,q-i)]';
        for j=1:size(Lambda,1),
          if lambdades(k) == Lambda(j,j),
            kident(k) = j; break
          else,
            kident(k) = 0;
          end
        end
        if kident(k) == 0
          G_OL_lambdades = evalfr(G_OL,lambdades(k));
          P(:,kP) = -inv(G_OL_lambdades)*e_vect;
          V(:,kP) = inv(A-lambdades(k)*eye(n))*B*P(:,kP);
        else
          disp('*** Warning : chosen eigenvalue = one open-loop eigenvalue ***')
          disp('*** Corresponding parameter vector = 0, and open-loop eigenvector kept ***'); beep
          P(:,kP) = 0; V(:,kP) = Vp(:,kident(k))
        end
        kP = kP+1;
        k=k+1;
        if i<q
            i=i+1;
        else
            i=1;
        end
      end
    end
    L_decoupling = P*inv(V);
  else
    disp(['Number of open-loop zeros is different frome n - Delta ==> '...
          'design impossible']); beep;
    L=0; M=0;
    return
  end
end  
L = Rounding_Matrix(L_decoupling); M = inv(Dstar)*diag([knum]);
