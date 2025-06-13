%% Convex Optimization

close all; clear; clc;

%% Gradient method

ALPHA = 0.01;   BETA = 0.5; MAXITERS = 1e+3;    GRADTOL = 1e-3;
n = 100; m = 200;

x = zeros(n,1);
A = randn(n, n); % rand(n,n)

k1 = 0; k2 = 0; kk1 = [];   kk2 = [];
vals = [];  tt = [];

for iter = 1:MAXITERS
    val = -sum(log(1-A*x)) - sum(log(1+x)) - sum(log(1-x));
    grad = A'*(1./(1-A*x)) - 1./(1+x) + 1./(1-x);
    
    if norm(grad) < GRADTOL
        break;
    end
    
    v = -grad;
    fprime = grad'*v;
    t = 1;
    
    while ((max(A*(x+t*v)) >= 1) | (max(abs(x+t*v)) >= 1))
        t = BETA*t;
        k1 = k1+1;
        tt = [tt, t];
    end
    
    while ( -sum(log(1-A*(x+t*v))) - sum(log(1-(x+t*v).^2)) > ...
    val + ALPHA*t*fprime )
        t = BETA*t;
        k2 = k2 + 1;
        tt = [tt, t];
    end
    
    x = x+t*v;
    
    kk1 = [kk1, k1];
    kk2 = [kk2, k2];
    vals = [vals, val];
    
end

%% plot (a)

figure(1)
semilogy((vals - vals(end)))
xlabel('iters')
ylabel('{f(x^{(k)}) - p^*}')
title('Plot of obejctive function')
grid on

figure(2)
stem(tt, ':diamondr', 'filled')
xlabel('iters')
ylabel('t')
grid on
title('step length')

%% Newton's method

n = 100; m = 200;
ALPHA = 0.01;   BETA = 0.5; MAXITERS = 1000;    NTTOL = 1e-8;
x = zeros(n,1); A = randn(n, n); % rand(n, n)

vals = [];
t_vals = [];

for iter = 1: MAXITERS
    val = -sum(log(1-A*x)) - sum(log(1+x)) - sum(log(1-x));
    d = 1./(1-A*x);
    grad = A'*d - 1./(1+x) + 1./(1-x);
    hess = A'*diag(d.^2)*A + diag(1./(1+x).^2 + 1./(1-x).^2);
    v = -hess\grad;
    fprime = grad'*v;
    
    if abs(fprime) < NTTOL
        break;
    end
    
    t = 1;
    t_vals =[t_vals, t];
    
    while ((max(A*(x+t*v)) >= 1) | (max(abs(x+t*v)) >= 1))
        t = BETA*t;
        t_vals = [t_vals, t];
    end
    
    while ( -sum(log(1-A*(x+t*v))) - sum(log(1-(x+t*v).^2)) > ...
    val + ALPHA*t*fprime )
        t = BETA*t;
        t_vals = [t_vals, t];
    end
    
    x = x+t*v;
    vals = [vals, val];
end

%% plot (b)

figure(3)
semilogy((vals - vals(end)))
xlabel('iters')
ylabel('{f(x^{(k)}) - p^*}')
title('Plot of obejctive function')
grid on

figure(4)
stem(t_vals, ':diamondr', 'filled')
xlabel('iters')
ylabel('t')
grid on
title('step length')

%%