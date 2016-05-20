function [Sect , Sectn , M] = Sector_fun(n , A)

% Example
% [A,B,C,D] = zp2ss(1,[-1+10i -1-10i -3 -3 -12],1);
% [Sect , Sectn , M] = Sector_fun([4 2 ; 3 5] , A);
% eigA = eig(M(:,:,1)\A*M(:,:,1));
% n : number of section which can be a matrix

count = numel(n);

for k = 1 : count
    condi = 1;
    Qk = A;
   while(condi > 1e-6) % obtain matrix sector
        Qk1 = Qk * (2*eye(size(Qk, 1)) + (n(k) - 2) * Qk^n(k)) / (eye(size(Qk, 1)) + (n(k) - 1) * Qk^n(k));
        condi = norm(Qk1 -Qk);
        Qk = Qk1;
    end
    Sect(:,:,k) = Qk; % page
    for q = 0:n(k)-1 % obtain q-th matrix sector
        tmp = 0;
        for i = 1:n(k)
            tmp= tmp + (Qk*exp(-1j*2*pi*q/n(k)))^(i-1)/n(k);
        end
        Sectn(:,:,q+1, k) = tmp; % box
    end
end

% obtain independent eigenvector
for k=1 : count
    tmpM=[];
    for i=1:n(k)
        numberOfVector = sum(logical(eig(Sectn(:,:,i,k))>0.5));
        j=0;
        tmp = [];
        rankOld = 0;
        while(numberOfVector)
            j=j+1;
            tmp(:,j) = Sectn(:,j,i,k);
            if(rankOld == rank(tmp))
                tmp(:,j)=[];
            end
            rankOld = rank(tmp);
            numberOfVector = numberOfVector - 1;
        end
        tmpM = [tmpM tmp];
    end
    M(:, :, k) = tmpM;
end



