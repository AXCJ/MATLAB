clear;clc;close all
M = 12; % length
N = 2 ^ M - 1; % period
idx = [1 4 6 12]; % feedback bits
u_bit = [2:5]'; % continuous 4 bits for filter
A = [0.5828 0.8948 1.19 1.0]; % filter coef.
B = [-0.5828 -0.8948  -1.19]; % filter coef. (right and left are reversed )

reg = randi([0 1],1,12); % get original registor
reg(find(~reg))=-1; % replace 0 to -1

u(1:3,:) = 0; % initial value
for j = 1:1:2^12*4 % 4 period
    Mreg(j+3,:) = reg; % save the reg
    u(j+3,:) = A*reg(u_bit)' + B*u(j : j+2); % filter
    
    % software xor
    % determining number of 1 is odd or even
    % the result is 1 if number of 1 is odd.
    shift = mod(size(find(reg(idx)==1),2), 2); % the shift bit
    if mod(size(find(reg(idx)==1),2), 2) == 0
        shift = -1;
    end
    % software xor --end
    reg=[shift reg(1:11)]; % shift it
end

StandardDeviation = std(u(N+1:2*N),1)
Mean = mean(u(N+1:2*N))

% compute autocorrelation r and correlation matrix R
% r(12) mean r(0); r(1) mean r(-11)
for k = -11:11
    r(k+12) = (u(N+1 : 2*N)' * u(N+1-k : 2*N-k))/N;
end
% make correlation matrix R
for i = 0 :1: 11
    R(i+1,:) = r(12-i:end-i);
end
% compute autocorrelation r and correlation matrix R --end

% plot
plot(u(N+1 : 4*N))
axis([-inf inf -inf inf])
xlabel('step')
line([N N],[min(u(N+1 : 2*N)), max(u(N+1 : 2*N))],'color','r')
line([2*N 2*N],[min(u(N+1 : 2*N)), max(u(N+1 : 2*N))],'color','r')
text(N,max(u(N+1 : 2*N))+0.1,'4095','color','r')
text(2*N,max(u(N+1 : 2*N))+0.1,'8190','color','r')

% uu=[1 2 3 4]';
% k=0
% for k = 0:1
%     r(k+1) = ((uu(1 : 4)'-mean(uu(1 : 4))) * (uu(1-k : 4-k)-mean(uu(1-k : 4-k))))/4
% end
% (std(uu(4+1 : 2*4)))




randi

