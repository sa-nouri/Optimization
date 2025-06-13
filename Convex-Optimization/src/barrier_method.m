%% Newton's Method for implementaion a barrier method

close all; clear; clc;

%% Parameter Initializing

m = 100; n = 500;
rand('seed',0);     randn('seed',0);

A = [rand(m-1, n); ones(1, n)];
x_0 = rand(n,1) + 0.1;
b = A * x_0;
c = randn(n, 1);

%% Xcript that generates data and tests the functions std_form_LP_acent

% analytic centering
[x_star, nu_star, lambda_hist] = lp_acent(A, b, c, x_0);
figure(1)
semilogy(lambda_hist, 'bo-')
xlabel('iterations')
ylabel('\lambda^2/2')
grid on

%% std_form_LP_barrier and solve the LP with barrier

[x_star, history, gap] = lp_barrier(A, b, c, x_0);
[xx, yy] = stairs(cumsum(history(1,:)),history(2,:));

figure(2)
semilogy(xx,yy, 'bo-');
xlabel('iters')
ylabel('gap')
grid on

% solve LP using cvx for comparison
cvx_begin
    variable x(n)
    minimize(c' * x)
    subject to
        A*x == b;
        x >= 0;
cvx_end

fprintf('\n\nOptimal value found by barrier method:\n');
p_star = c'*x_star
fprintf('Optimal value found by CVX:\n');
cvx_optval
fprintf('Duality gap from barrier method:\n');
gap

%%  We test our LP solver on two problem instances, one infeasible, and one feasible
% solves standard form LP for two problem instances

m = 100;    n = 500;

% infeasible problem instance
rand('state', 0);   randn('state', 0);
A = [randn(m-1, n); ones(1, n)];
b = randn(m, 1); c = randn(n, 1);

[x_star, p_star, gap, status, nsteps] = lp_solve(A, b, c);

% solve LP using cvx for comparison
cvx_begin
    variable x(n)
    minimize(c' * x)
    subject to
        A*x == b;
        x >= 0;
cvx_end


% feasible problem instance
A = [randn(m-1, n); ones(1, n)];
v = rand(n, 1) + 0.1;
b = A*v;
c = randn(n, 1);
[x_star, p_star, gap, status, nsteps] = lp_solve(A, b, c);

% solve LP using cvx for comparison
cvx_begin
    variable x(n)
    minimize(c'*x)
    subject to
        A*x == b;
        x >= 0;
cvx_end

fprintf('\n\nOptimal value found by barrier method:\n');
p_star

fprintf('Optimal value found by CVX:\n');
cvx_optval

fprintf('Duality gap from barrier method:\n');
gap

%% Defining the lp acent function
% The following functions compute the analytic center using Newton’s method
function [x_star, nu_star, lambda_hist] = lp_acent(A,b,c,x_0)
    % solves problem
    % minimize c’*x - sum(log(x))
    % subject to A*x = b
    % using Newton’s method, given strictly feasible starting point x0
    % input (A, b, c, x_0)
    % returns primal and dual optimal points
    % lambda_hist is a vector showing lambda^2/2 for each newton step
    % returns [], [] if MAXITERS reached, or x_0 not feasible

    % algorithm parameters
    ALPHA = 0.01;   BETA = 0.5;     EPSILON = 1e-3; MAXITERS = 100;
    
    if (min(x_0) <= 0) || (norm(A*x_0 - b) > 1e-3) % x0 not feasible
        fprintf('FAILED');
        nu_star = []; x_star = []; lambda_hist=[];
        return;
    end
    
    m = length(b);  n = length(x_0);    x = x_0; lambda_hist = [];
    
    for iter = 1 : MAXITERS
        H = diag(x.^(-2));
        g = c - x.^(-1);
        % lines below compute newton step via whole KKT system
        % M = [ H A’; A zeros(m,m)];
        % d = M\[-g; zeros(m,1)];
        % dx = d(1:n);
        % w = d(n+1:end);
        % newton step by elimination method
        
        w = (A*diag(x.^2)*A')\(-A*diag(x.^2)*g);
        dx = -diag(x.^2)*(A'*w + g);
        lambdasqr = -g'*dx; % dx’*H*dx;
        lambda_hist = [lambda_hist lambdasqr/2];
        if lambdasqr/2 <= EPSILON 
            break;
        end
        
        % backtracking line search
        % first bring the point inside the domain
        
        t = 1;
        while min(x+t*dx) <= 0
            t = BETA*t;
        end
        
        % now do backtracking line search
        while c'*(t*dx)-sum(log(x+t*dx))+sum(log(x))-ALPHA*t*g'*dx> 0
            t = BETA*t;
        end
        x = x + t*dx;
    end
    if iter == MAXITERS % MAXITERS reached
    fprintf('ERROR: MAXITERS reached.\n');
    x_star = []; nu_star = [];
    else
    x_star = x;
    nu_star = w;
    end
end

%% Defining the lp barrier function
% The following functions solve the LP using the barrier method
function [x_star, history, gap] = lp_barrier(A, b, c, x_0)
    % solves standard form LP
    % minimize c^T x
    % subject to Ax = b, x >=0;
    % using barrier method, given strictly feasible x0
    % uses function std_form_LP_acent() to carry out centering steps
    % returns:
    % - primal optimal point x_star
    % - history, a 2xk matrix that returns number of newton steps
    % in each centering step (top row) and duality gap (bottom row)
    % (k is total number of centering steps)
    % - gap, optimal duality gap
    
    % barrier method parameters
    T_0 = 1;    MU = 20;
    EPSILON = 1e-3; % duality gap stopping criterion
    n = length(x_0);    t = T_0;    x = x_0;
    
    history = [];
    while(1)
        gap = n/t;
        if gap < EPSILON
            break;
        end
        [x_star, nu_star, lambda_hist] = lp_acent(A, b, t*c, x);
        x = x_star;
        
        history = [history [length(lambda_hist); gap]];
        t = MU*t;
    end
end

%% Defining lp solving
% The following functions implement the full LP solver (phase I and phase II)

function [x_star, p_star, gap, status, nsteps] = lp_solve(A, b, c)
    % solves the LP
    % minimize c^T x
    % subject to Ax = b, x >= 0;
    % using a barrier method
    % computes a strictly feasible point by carrying out
    % a phase I method
    % returns:
    % - a primal optimal point x_star
    % - the primal optimal value p_star
    % - status: either ’Infeasible’ or ’Solved’
    % - nsteps(1): number of newton steps for phase I
    % - nsteps(2): number of newton steps for phase II
    [m,n] = size(A);    nsteps = zeros(2, 1);
    
    % phase I
    x0 = A\b; t0 = 2 + max(0, -min(x0));
    A1 = [A, -A * ones(n, 1)];
    b1 = b - A * ones(n, 1);
    z0 = x0 +t0 * ones(n, 1) - ones(n, 1);
    c1 = [zeros(n, 1); 1];
    [z_star, history, gap] = lp_barrier(A1, b1, c1, [z0; t0]);
    if (z_star(n+1) >= 1)
        fprintf('\nProblem is infeasible\n');
        x_star = []; p_star = Inf; status = 'Infeasible';
        nsteps(1) = sum(history(1, :)); gap = [];
        return;
    end
    
    fprintf('\nFeasible point found\n');
    nsteps(1) = sum(history(1, :));
    x_0 = z_star(1:n) - z_star(n+1) * ones(n, 1) + ones(n, 1);
    
    % phase II
    [x_star, history, gap] = lp_barrier(A, b, c, x_0);
    status = 'Solved';
    p_star = c'*x_star;
    nsteps(2) = sum(history(1, :));
    
end

