%% Brute force solution for solving a subset of convex constraints

close all; clear; clc;

%%

% Run the given file for loading the value of variables.
run satisfy_some_constraints_data.m;

%% Simulation

feasible_threshold = 1e-5;

cvx_begin
    variables x(n) mu_opt
    minimize(c' * x)
    subject to
        sum(pos(mu_opt + A * x - b)) <= mu_opt * (m - k);
        mu_opt >= 0;
cvx_end

fprintf('The optimal value of lambda is: %f\n', 1 / mu_opt)
fprintf('The objective value: %f\n', cvx_optval)
fprintf('The number of satisfied constraints: %d\n', nnz(A * x - b <= feasible_threshold))

[~, Indexes] = sort(A * x - b);
selected_constraints = Indexes(1:k);

cvx_begin
    variables x(n)
    minimize(c' * x)
    subject to
        A(selected_constraints, :) * x <= b(selected_constraints);
cvx_end

fprintf('The objective value for selected constraints: %f\n', cvx_optval)

%%