%% Risk-return trade-off in portfolio optimization

close all; clear; clc;

%% 

n=4;
novals = 100;   p_bar = [.12 .10 .07 .03]';
Sigma = [0.0064 0.0008 -0.0011 0;
         0.0008 0.0025 0 0;
        -0.0011 0 0.0004 0;
         0 0 0 0];

returns = zeros(1, novals);  std_values = zeros(1, novals);
pfs = zeros(n, novals);      mu_values = logspace(0, 7, novals);

i = 1;

while i <= novals
    mu = mu_values(i);
    
    cvx_begin
        variable x(n);
        minimize (-p_bar' * x + mu * quad_form(x, Sigma));
        subject to
            sum(x) == 1;
            x >= 0;
    cvx_end
    
    returns(i) = p_bar' * x;
    std_values(i) = sqrt(x' * Sigma * x);
    pfs(1:n, i)= x;
    i = i +1;
end

figure(1)
plot(std_values, returns)
grid on
grid minor
xlabel('Standard Deviation')
ylabel('Returns')

figure(2)
plot(std_values, pfs(1,:)')
hold on
plot(std_values, (pfs(1,:)+ pfs(2, :))')
plot(std_values, (pfs(1,:) + pfs(2,:)+ pfs(3,:))')
grid on
grid minor
xlabel('Standard Deviation')
ylabel('Pfs')
hold off

%% Apply loss Constraint to portfolio optimization

betha = 0; novals = 100;
etas = logspace(-4, -1, novals);

[V, D] = eig(Sigma);
sqrt_Sigma = V*diag(sqrt(diag(D))) * V';

pfs = zeros(n, novals);     returns = zeros(1, novals);

k = 1;
while k <= novals
    eta = etas(k);
    gamma = -(2^1/2) * erfcinv(2 * eta);
    
    cvx_begin
        variable x(n)
        maximize(p_bar' * x)
        subject to
            x >= 0;
            sum(x) == 1;
            p_bar' * x + gamma * norm(sqrt_Sigma * x, 2) >= betha;
    cvx_end
    
    returns(k) = p_bar' * x;
    pfs(:, k) = x;
    k = k +1;
end


figure(1)
semilogx(etas, returns)
xlabel('Eta')
ylabel('Returns')
grid on
grid minor

figure(2)
semilogx(etas, pfs(1,:)') 
hold on
semilogx(etas, (pfs(1, :) + pfs(2, :))')
semilogx(etas, (pfs(1, :) + pfs(2, :) + pfs(3, :))')
xlabel('Eta')
ylabel('Pfs')
grid on
grid minor
hold off

%% Monte carlo simulation

N=10000;    eta = 0.05;
randn('state', 0);
gamma = sqrt(2) * erfcinv(2 * (1 - eta));
returns = [];

cvx_begin
    variable x(n)
    maximize(p_bar' * x)
    subject to
        p_bar'*x + gamma * norm(sqrt_Sigma * x, 2) >= 0;
        sum(x) == 1;
        x >= 0;
cvx_end


i = 1;
while i < N
    p = p_bar + sqrt_Sigma * randn(n, 1);
    returns = [returns p'*x];
    i = i + 1;
end


figure(1)
histogram(returns, 50)
hold on
grid on
grid minor
plot([0;0],[0;700],'r--')
title('Histogram of Returns')

disp('The emprical mean of returns \n')
disp(mean(returns))

disp('Percent of positive returns are \n')
disp(sum(returns < 0)/N)
%%