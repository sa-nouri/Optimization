%% Generator Dispatch

close all; clear; clc;

%% solve generator dispatch problem
run gen_dispatch_data

cvx_begin quiet
    variable p(n, T)
    dual variable Q
    p <= Pmax' * ones(1, T);
    p >= Pmin' * ones(1, T);
    Q: sum(p, 1) == d;
    abs(p(:, 2:T) - p(:, 1:T-1)) <= R' * ones(1, T-1);
    power_cost = sum(alpha * p + beta *(p.^2));
    change_power_cost = sum(gamma * abs(p(:, 2:T) - p(:, 1:T-1)));
    minimize (power_cost+change_power_cost)
cvx_end

subplot(3, 1, 1)
plot(t, d)
title('demand')
grid on
grid minor

subplot(3, 1, 2)
plot(t, p);
title('generator powers')
grid on
grid minor

subplot(3, 1, 3)
plot(t, Q);
title('power prices')
grid on
grid minor

%% Solve problem for each generator

cvx_begin quiet
    variable pp(n, T)
    pp <= Pmax' * ones(1, T);
    pp >= Pmin' * ones(1, T);
    %sum(pp,1) == d;
    abs(pp(:, 2:T) - pp(:, 1:T-1)) <= R' * ones(1, T-1);
    power_cost = sum(alpha * pp + beta *(pp.^2));
    change_power_cost = sum(gamma*abs(pp(:,2:T)-pp(:,1:T-1)));
    minimize (power_cost+change_power_cost-sum(pp*Q'))
cvx_end

rel_error = norm(pp-p)/norm(p)