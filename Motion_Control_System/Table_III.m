
clear;clc;close all;
ht = zeros(8,3);
for j = 1:3
N = input('Please enter N:');
M = input('Please enter M:');
    A = [];
    qt = [];


    for i = 1 : M
        A = [A ;i];
    end
    columnA = A;
    for i = 2 : N
        A = [A columnA.^i];
    end
    A = [ones(size(A,1),1) A]
    for i = 0 : N
        qt = [qt i*M^(i-1)];
    end
    qt
    ht(1:size((qt*pinv(A))',1),j) = (qt*pinv(A))';
end
ht
colnames = {'LSF 1/4', 'LSF 2/8', 'LSF 3/8'};
rows = {'h1','h2','h3','h4','h5','h6','h7','h8'};
t = uitable('Data',ht,'ColumnName',colnames,'RowName',rows,'ColumnWidth',{65},'FontSize',12);