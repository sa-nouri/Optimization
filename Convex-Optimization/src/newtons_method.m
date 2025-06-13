%% Other versions of Newton's method implementation

close all; clear; clc;

%% Re-using the Hessian

m = 200; n = 100;
rand('seed', 0);     randn('seed', 0);

ALPHA = 0.01;   BETA = 0.5; MAXITERS = 1000;    NTTOL = 1e-8;
x = zeros(n,1); A = rand(n, n);
N = 8;

vals = [];
t_vals = [];

for iter = 0: MAXITERS
    val = -sum(log(1-A*x)) - sum(log(1+x)) - sum(log(1-x));
    d = 1./(1-A*x);
    grad = A'*d - 1./(1+x) + 1./(1-x);
    if (mod(iter, N) == 0)
        hess = A' * diag(d.^2) *A + diag(1./(1+x).^2 + 1./(1-x).^2);
    end
    v = -hess\grad;
    fprime = grad' * v;
    
    if abs(fprime) < NTTOL
        break;
    end
    
    t = 1;
    
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

%% Saving values

% xr = []
% xr = vals; N =1
% xr1 = vals; N= 2
% xr2 = vals; N = 5;
% xr3 = vals; N =8;

%% plot(a)

figure(1)
semilogy((xr - xr(end)))
xlabel('iters')
ylabel('{f(x^{(k)}) - p^*}')
title('Plot of obejctive function for different values of N')
grid on
hold on
semilogy((xr - xr(end)))
semilogy((xr1- xr1(end)))
semilogy((xr2 - xr2(end)))
semilogy((xr3 - xr3(end)))
legend('N=1', 'N=2', 'N=5', 'N=8')
hold off

%% Diagonal approximation

vals = [];
t_vals = [];

for iter = 1: MAXITERS
    val = -sum(log(1-A*x)) - sum(log(1+x)) - sum(log(1-x));
    d = 1./(1-A*x);
    grad = A'*d - 1./(1+x) + 1./(1-x);
%     hess = A'*diag(d.^2)*A + diag(1./(1+x).^2 + 1./(1-x).^2);
    hess = diag(1./(1+x).^2 + 1./(1-x).^2);
    v = -hess\grad;
    fprime = grad'*v;
    
    if abs(fprime) < NTTOL
        break;
    end
    
    t = 1;
    
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

%% plot(b)

figure(2)
semilogy((vals - vals(end)))
xlabel('iters')
ylabel('{f(x^{(k)}) - p^*}')
title('Plot of obejctive function')
grid on

