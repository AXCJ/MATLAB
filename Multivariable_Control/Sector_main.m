clear; clc;

[A,B,C,D] = zp2ss(1,[-1+9i -1-9i -2 -11 10],1);
[Sect1 , Sectn1 , M1] = Sector_fun([4] , A);
for i = 1 : size(M1 , 3)
%     eig(M(:,:,i)\A*M(:,:,i))
    Astar1 = M1(:,:,i)\A*M1(:,:,i);
    eig(Astar1);
end
[Sect2 , Sectn2 , M2] = Sector_fun([2] , A+10*eye(size(A,1)));
for i = 1 : size(M2 , 3)
    Astar2 = M2(:,:,i)\A*M2(:,:,i);
    eig(Astar2);
end

Astar(:,:,2) = A;
Astar(:,:,1) = A+10*eye(size(A,1));
[Sect3 , Sectn3 , M3] = Sector_fun1([2 4] , Astar , [1 3]);
for i = 1 : size(M3 , 3)
    Astar3 = M3(:,:,i)\A*M3(:,:,i);
    eig(Astar3)
end
